# Drug Concentration in Blood

A lightweight Flask web app that simulates and visualizes drug concentration in blood plasma over time using a continuous-time Linear Time-Invariant (LTI) pharmacokinetic model.

## Features

- Simulates oral dosing and continuous IV infusion curves
- Supports preset drugs and custom parameter entry
- Displays therapeutic window, peak concentration, and safety status
- Includes interpretation messages for toxic/subtherapeutic/valid dosing

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

## Notes

- Matplotlib is configured to use a non-GUI backend (`Agg`) so plotting works safely in Flask request threads.
- This is a development server setup; use a production WSGI server for deployment.
