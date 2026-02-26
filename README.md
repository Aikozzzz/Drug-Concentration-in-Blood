# Drug Concentration in Blood

A lightweight Flask web app that simulates and visualizes drug concentration in blood plasma over time using a continuous-time Linear Time-Invariant (LTI) pharmacokinetic model.

## Features

- Simulates oral dosing and continuous IV infusion curves
- Interactive hover-enabled graph with per-point explanation (below effective / therapeutic / toxic)
- Supports preset drugs and custom parameter entry
- Displays therapeutic window, peak concentration, and safety status
- Includes interpretation messages for toxic/subtherapeutic/valid dosing
- Includes a built-in virtual assistant chat for explaining graph behavior and parameter effects
- Chat opens from a top-right icon and can be closed again
- Rule-based assistant uses predefined question choices (no free-text typing)

## Project Structure

- `app.py` - Flask backend and simulation logic
- `templates/index.html` - UI form and graph display
- `static/` - static assets

## Requirements

- Python 3.10+
- pip

## Installation

```bash
python -m venv .venv
```

### Windows (PowerShell)

```powershell
.\.venv\Scripts\activate
pip install -r requirements.txt
```

### macOS/Linux

```bash
source .venv/bin/activate
pip install -r requirements.txt
```

## Run the App

```bash
python app.py
```

Open in browser:

`http://127.0.0.1:5000`

## Virtual Assistant Chat

- Click the top-right chat icon (💬) to open the assistant panel.
- Choose a predefined question and click **Ask**.
- Supported topics include status, dose adjustment, interval effect, `ka`, `ke`, `Vd`, and oral vs IV interpretation.
- The assistant uses current simulation context (including status) for responses.

## Interactive Graph

- The main pharmacokinetic chart is interactive.
- Hover over oral/IV curve points to see time, concentration, and zone explanation.
- Threshold lines for toxic and minimum effective levels are also hoverable.

## Notes

- Matplotlib is configured to use a non-GUI backend (`Agg`) so plotting works safely in Flask request threads.
- Plotly is loaded in the browser via CDN for interactive chart rendering.
- This is a development server setup; use a production WSGI server for deployment.
