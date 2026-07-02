from flask import Flask, render_template, request
from optimizer import optimize_both

app = Flask(__name__)


@app.route("/")
def index():
    return render_template(
        "index.html",
        default_weights=None,
        seq_type="dna",
        matrix_name=None,
    )


@app.route("/align", methods=["POST"])
def align():
    seq1 = request.form.get("seq1", "").strip()
    seq2 = request.form.get("seq2", "").strip()

    file1 = request.files.get("file1")
    file2 = request.files.get("file2")

    if file1 and file1.filename:
        raw = file1.read().decode("utf-8")
        lines = raw.splitlines()
        if lines and lines[0].startswith(">"):
            lines = lines[1:]
        seq1 = "".join(lines).strip()
    if file2 and file2.filename:
        raw = file2.read().decode("utf-8")
        lines = raw.splitlines()
        if lines and lines[0].startswith(">"):
            lines = lines[1:]
        seq2 = "".join(lines).strip()

    if not seq1 or not seq2:
        return render_template("index.html", error="Please provide both sequences.")

    seq_type = request.form.get("seq_type", "dna")
    matrix_name = request.form.get("matrix_name") or None

    try:
        match_w = (
            float(request.form["match_weight"])
            if request.form.get("match_weight") and seq_type == "dna"
            else None
        )
        mismatch_w = (
            float(request.form["mismatch_weight"])
            if request.form.get("mismatch_weight") and seq_type == "dna"
            else None
        )
        gap_w = (
            float(request.form["gap_weight"])
            if request.form.get("gap_weight")
            else None
        )
    except ValueError:
        return render_template("index.html", error="Invalid weight value.")

    result = optimize_both(
        seq1,
        seq2,
        match_w,
        mismatch_w,
        gap_w,
        seq_type=seq_type,
        matrix_name=matrix_name,
    )
    default_weights = result["needleman_wunsch"]["weights"]
    return render_template(
        "index.html",
        result=result,
        seq1=seq1,
        seq2=seq2,
        seq1_preview=seq1[:60],
        seq2_preview=seq2[:60],
        default_weights=default_weights,
        seq_type=seq_type,
        matrix_name=matrix_name,
    )


if __name__ == "__main__":
    app.run(debug=True)
