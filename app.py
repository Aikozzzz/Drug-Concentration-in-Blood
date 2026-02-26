from flask import Flask, render_template, request, jsonify
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import io
import base64
import json
import re
from scipy import signal

app = Flask(__name__)

# ==============================
# Updated Drug Database
# ==============================

DRUG_DB = {
    "Theophylline": {
        "ka": 1.0,
        "ke": 0.1,
        "dose": 200,
        "doses_per_day": 4,
        "interval": 16 / 4,
        "toxic": 20,
        "min_eff": 5,
        "Vd": 35,
        "F": 1.0,
        "scale": 1
    },
    "Antibiotics": {   # generic beta-lactam profile
        "ka": 1.3,
        "ke": 0.5,
        "dose": 500,
        "doses_per_day": 3,
        "interval": 16 / 3,
        "toxic": 50,
        "min_eff": 10,
        "Vd": 20,
        "F": 0.9,
        "scale": 1
    },
    "Furosemide": {
        "ka": 1.2,
        "ke": 0.35,
        "dose": 40,
        "doses_per_day": 2,
        "interval": 16 / 2,
        "toxic": 8,
        "min_eff": 1,
        "Vd": 12,
        "F": 0.7,
        "scale": 1
    },
    "Metronidazole": {
        "ka": 1.1,
        "ke": 0.15,
        "dose": 400,
        "doses_per_day": 3,
        "interval": 16 / 3,
        "toxic": 40,
        "min_eff": 8,
        "Vd": 40,
        "F": 1.0,
        "scale": 1
    },
    "Ciprofloxacin": {
        "ka": 1.5,
        "ke": 0.25,
        "dose": 500,
        "doses_per_day": 3,
        "interval": 16 / 3,
        "toxic": 6,
        "min_eff": 1,
        "Vd": 30,
        "F": 0.8,
        "scale": 1
    },
    "Custom": {
        "ka": 1.2,
        "ke": 0.3,
        "dose": 100,
        "doses_per_day": 3,
        "interval": 16 / 3,
        "toxic": 50,
        "min_eff": 5,
        "Vd": 40,
        "F": 1.0,
        "scale": 1
    }
}

SIM_TIME = 16  # 16 hours simulation


def parse_float(form, key, default):
    try:
        return float(form.get(key, default))
    except:
        return float(default)


def validate_params(params):
    errors = []

    if params["ka"] <= 0:
        errors.append("ka must be greater than 0.")
    if params["ke"] <= 0:
        errors.append("ke must be greater than 0.")
    if params["Vd"] <= 0:
        errors.append("Vd must be greater than 0.")
    if params["dose"] <= 0:
        errors.append("Dose must be greater than 0.")
    if params["doses_per_day"] < 1:
        errors.append("Doses/day must be at least 1.")
    if params["interval"] <= 0:
        errors.append("Interval must be greater than 0.")
    if params["toxic"] <= 0:
        errors.append("Toxic limit must be greater than 0.")
    if params["min_eff"] < 0:
        errors.append("Minimum effective concentration cannot be negative.")
    if params["toxic"] <= params["min_eff"]:
        errors.append("Toxic limit must be greater than minimum effective concentration.")
    if params["delay2"] < 0 or params["delay3"] < 0:
        errors.append("Dose delays cannot be negative.")

    return errors


def build_assistant_reply(message, context, intent=None):
    text = (message or "").strip().lower()
    status = str(context.get("status") or "").upper()
    tokens = set(re.findall(r"[a-zA-Z0-9_]+", text))

    if intent == "status":
        if status == "TOXIC":
            return "Your current status is TOXIC. Consider reducing dose, increasing interval, or both to lower peak concentration."
        if status == "SUBTHERAPEUTIC":
            return "Your current status is SUBTHERAPEUTIC. Consider increasing dose or shortening interval to reach effective levels."
        if status == "SAFE":
            return "Your current status is SAFE based on current thresholds. You can still check trough behavior if dose timing changes."
        if status == "INVALID INPUT":
            return "Inputs are invalid. Please fix highlighted issues first (for example non-positive dose, interval, or Vd)."
        return "Run a simulation first, then ask about status and I will interpret it for you."

    if intent == "dose_adjust":
        return "Higher dose or more frequent dosing increases concentration and peak risk. Lowering either usually reduces toxicity risk."

    if intent == "interval_effect":
        return "Longer intervals usually reduce accumulation but can cause subtherapeutic dips. Shorter intervals increase stability but may raise toxicity risk."

    if intent == "ka_meaning":
        return "Higher ka means faster absorption and earlier/higher peaks. Lower ka smooths and delays peak concentration."

    if intent == "ke_meaning":
        return "Higher ke means faster elimination, lower accumulation, and lower troughs. Lower ke increases drug persistence."

    if intent == "vd_meaning":
        return "Larger Vd generally lowers plasma concentration for the same dose, while smaller Vd raises concentration."

    if intent == "oral_vs_iv":
        return "The dashed IV curve represents continuous infusion, while the oral curve shows pulse-like dose events with absorption and elimination dynamics."

    if not text:
        return "Please type a question about dosing, toxicity, or the graph."

    if text in {"hello", "hi", "hey"}:
        return "Hi! I can help explain your graph, status, dose interval, and therapeutic window."

    if "status" in text or "safe" in text or "toxic" in text or "subtherapeutic" in text:
        if status == "TOXIC":
            return "Your current status is TOXIC. Consider reducing dose, increasing interval, or both to lower peak concentration."
        if status == "SUBTHERAPEUTIC":
            return "Your current status is SUBTHERAPEUTIC. Consider increasing dose or shortening interval to reach effective levels."
        if status == "SAFE":
            return "Your current status is SAFE based on current thresholds. You can still check trough behavior if dose timing changes."
        if status == "INVALID INPUT":
            return "Inputs are invalid. Please fix highlighted issues first (for example non-positive dose, interval, or Vd)."
        return "Run a simulation first, then ask about status and I will interpret it for you."

    if "interval" in text or "timing" in text or "delay" in text:
        return "Longer intervals usually reduce accumulation but can cause subtherapeutic dips. Shorter intervals increase stability but may raise toxicity risk."

    if "dose" in text or "doses/day" in text or "doses per day" in text:
        return "Higher dose or more frequent dosing increases concentration and peak risk. Lowering either usually reduces toxicity risk."

    if "ka" in tokens or "absorption" in tokens:
        return "Higher ka means faster absorption and earlier/higher peaks. Lower ka smooths and delays peak concentration."

    if "ke" in tokens or "elimination" in tokens:
        return "Higher ke means faster elimination, lower accumulation, and lower troughs. Lower ke increases drug persistence."

    if "vd" in tokens:
        return "Larger Vd generally lowers plasma concentration for the same dose, while smaller Vd raises concentration."

    if "iv" in text:
        return "The dashed IV curve represents continuous infusion, while the oral curve shows pulse-like dose events with absorption and elimination dynamics."

    return "I can explain status, dose interval, ka/ke/Vd effects, toxicity risk, and how to adjust schedule parameters."


def generate_plot(params):
    validation_errors = validate_params(params)
    if validation_errors:
        raise ValueError("; ".join(validation_errors))

    t = np.linspace(0, SIM_TIME, 1000)
    dt = t[1] - t[0]
    scale = params.get("scale", 1)

    # Transfer Functions
    num_oral = [params["F"] * params["ka"] / params["Vd"]]
    den_oral = [1, params["ka"] + params["ke"], params["ka"] * params["ke"]]
    sys_oral = signal.TransferFunction(num_oral, den_oral)

    num_iv = [1 / params["Vd"]]
    den_iv = [1, params["ke"]]
    sys_iv = signal.TransferFunction(num_iv, den_iv)

    # ========================
    # ORAL DOSE TIMING
    # ========================

    doses_per_day = int(params["doses_per_day"])
    interval = params["interval"]

    dose_times = []
    for i in range(doses_per_day):
        time_point = i * interval
        if i == 1:
            time_point += params["delay2"]
        elif i == 2:
            time_point += params["delay3"]
        dose_times.append(time_point)

    # Impulse train
    u_pill = np.zeros_like(t)
    for time_point in dose_times:
        if 0 <= time_point <= SIM_TIME:
            idx = np.argmin(np.abs(t - time_point))
            u_pill[idx] = params["dose"] / dt

    _, y_pill, _ = signal.lsim(sys_oral, U=u_pill, T=t)

    # ========================
    # IV INFUSION (constant)
    # ========================
    u_iv = np.full_like(t, params["dose"] / interval)
    _, y_iv, _ = signal.lsim(sys_iv, U=u_iv, T=t)

    y_pill *= scale
    y_iv *= scale

    toxic_scaled = params["toxic"] * scale
    min_eff_scaled = params["min_eff"] * scale

    # Plot
    plt.figure(figsize=(10, 5))
    plt.plot(t, y_pill, label="Oral doses", linewidth=2)
    plt.plot(t, y_iv, linestyle="--", label="Continuous IV")
    plt.axhline(y=toxic_scaled, linestyle=":", label="Toxic")
    plt.axhline(y=min_eff_scaled, linestyle=":", label="Min Effective")
    plt.fill_between(t, min_eff_scaled, toxic_scaled, alpha=0.1)

    plt.title(f"Pharmacokinetics - {params['name']}")
    plt.xlabel("Time (hours)")
    plt.ylabel("Plasma Concentration (mg/L)")
    plt.legend()
    plt.grid(alpha=0.2)

    img = io.BytesIO()
    plt.savefig(img, format="png", bbox_inches="tight")
    img.seek(0)
    plot_url = base64.b64encode(img.getvalue()).decode()
    plt.close()

    max_conc = float(np.max(y_pill))
    min_conc = float(np.min(y_pill))

    if max_conc > toxic_scaled:
        status = "TOXIC"
    elif max_conc < min_eff_scaled:
        status = "SUBTHERAPEUTIC"
    else:
        status = "SAFE"
        
     # =============================
    # USER-FRIENDLY MESSAGES
    # =============================
    messages = []

    # Toxicity message
    if max_conc > toxic_scaled:
        messages.append(
            "⚠ Drug concentration exceeds the toxic threshold. This schedule may cause adverse effects."
        )

    # Under-dose message
    if max_conc < min_eff_scaled:
        messages.append(
            "⚠ Drug levels never reach the minimum effective concentration. Therapeutic effect may not occur."
        )

    # Drops below effective level before next dose
    # Only warn here when the overall status is not SAFE
    if status != "SAFE" and min_conc < min_eff_scaled and max_conc > min_eff_scaled:
        messages.append(
            "⚠ Drug concentration falls below effective level before the next dose. Consider reducing interval."
        )

    # Delay warnings
    if params["delay2"] > 0.5:
        messages.append(
            "⏳ There is a noticeable delay in the second dose. This may reduce overall drug stability."
        )

    if params["delay3"] > 0.5:
        messages.append(
            "⏳ There is a delay in the third dose. Irregular dosing can cause concentration fluctuations."
        )

    # Too frequent dosing
    if interval < (SIM_TIME / doses_per_day) * 0.7:
        messages.append(
            "⚠ Doses are administered very frequently. Monitor for drug accumulation."
        )

    if not messages:
        messages.append(
            "✔ Drug concentration remains within the therapeutic window. Current dosing appears appropriate."
        )

    def point_explanation(concentration):
        if concentration > toxic_scaled:
            return "Above toxic threshold"
        if concentration < min_eff_scaled:
            return "Below minimum effective level"
        return "Within therapeutic window"

    oral_explanations = [point_explanation(float(v)) for v in y_pill]
    iv_explanations = [point_explanation(float(v)) for v in y_iv]

    plot_payload = {
        "name": params["name"],
        "t": [float(v) for v in t],
        "y_oral": [float(v) for v in y_pill],
        "y_iv": [float(v) for v in y_iv],
        "toxic": float(toxic_scaled),
        "min_eff": float(min_eff_scaled),
        "oral_explanations": oral_explanations,
        "iv_explanations": iv_explanations,
    }

    return plot_url, status, round(max_conc, 2), messages, plot_payload


@app.route("/", methods=["GET", "POST"])
def index():
    selected_drug = "Theophylline"
    base = DRUG_DB[selected_drug]

    form_values = {
        "ka": base["ka"],
        "ke": base["ke"],
        "Vd": base["Vd"],
        "dose": base["dose"],
        "doses_per_day": base["doses_per_day"],
        "interval": base.get("interval", 16 / base["doses_per_day"]),
        "delay2": 0,
        "delay3": 0,
        "toxic": base["toxic"],
        "min_eff": base["min_eff"],
    }

    plot_url = None
    plot_payload = None
    status = None
    max_c = None
    messages = None

    if request.method == "POST":
        selected_drug = request.form.get("drug_select", "Theophylline")
        if selected_drug not in DRUG_DB:
            selected_drug = "Theophylline"
        base = DRUG_DB[selected_drug]

        p = {
            "ka": parse_float(request.form, "ka", base["ka"]),
            "ke": parse_float(request.form, "ke", base["ke"]),
            "dose": parse_float(request.form, "dose", base["dose"]),
            "doses_per_day": int(parse_float(request.form, "doses_per_day", base["doses_per_day"])),
            "interval": parse_float(request.form, "interval", base.get("interval", 16 / base["doses_per_day"])),
            "delay2": parse_float(request.form, "delay2", 0),
            "delay3": parse_float(request.form, "delay3", 0),
            "toxic": parse_float(request.form, "toxic", base["toxic"]),
            "min_eff": parse_float(request.form, "min_eff", base["min_eff"]),
            "Vd": parse_float(request.form, "Vd", base["Vd"]),
            "F": base["F"],
            "scale": base["scale"],
            "name": selected_drug,
        }

        form_values = p
        input_errors = validate_params(p)
        if input_errors:
            status = "INVALID INPUT"
            messages = [f"⚠ {error}" for error in input_errors]
            plot_url = None
            plot_payload = None
            max_c = None
        else:
            plot_url, status, max_c, messages, plot_payload = generate_plot(p)

    return render_template(
        "index.html",
        plot_url=plot_url,
        plot_payload=plot_payload,
        status=status,
        max_c=max_c,
        messages=messages,
        drugs=list(DRUG_DB.keys()),
        selected=selected_drug,
        form_values=form_values,
        drugs_json=json.dumps(DRUG_DB)
    )


@app.route("/assistant-chat", methods=["POST"])
def assistant_chat():
    data = request.get_json(silent=True) or {}
    message = str(data.get("message", ""))
    intent = str(data.get("intent", "")).strip().lower() or None
    context = data.get("context", {})
    if not isinstance(context, dict):
        context = {}

    reply = build_assistant_reply(message, context, intent=intent)
    return jsonify({"reply": reply})


if __name__ == "__main__":
    app.run(debug=True)
