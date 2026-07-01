from flask import Flask, render_template, request
from optimizer import optimize_both

app = Flask(__name__)


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/align", methods=["POST"])
def align():
    seq1 = request.form.get("seq1", "").strip()
    seq2 = request.form.get("seq2", "").strip()

    file1 = request.files.get("file1")
    file2 = request.files.get("file2")

    if file1 and file1.filename:
        seq1 = file1.read().decode("utf-8").replace("\n", "").replace("\r", "").strip()
    if file2 and file2.filename:
        seq2 = file2.read().decode("utf-8").replace("\n", "").replace("\r", "").strip()

    if not seq1 or not seq2:
        return render_template("index.html", error="Please provide both sequences.")

    result = optimize_both(seq1, seq2)
    return render_template(
        "index.html", result=result, seq1_preview=seq1[:60], seq2_preview=seq2[:60]
    )


if __name__ == "__main__":
    app.run(debug=True)
