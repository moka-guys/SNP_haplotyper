from argparse import Namespace
from datetime import datetime
import excel_parser
from flask import Flask, render_template, Response, send_file, jsonify, request, session
from flask_wtf import FlaskForm
import os
import time
import random
from wtforms import FileField, SubmitField, MultipleFileField
from werkzeug.utils import secure_filename


app = Flask(__name__)
app.config["UPLOAD_FOLDER"] = "./uploads/tmp_folder"
app.config["SECRET_KEY"] = "catchmeifyoucan"
app.config["UPLOAD_EXTENSIONS"] = [".txt", ".csv", ".xlsm", ".xlsx"]
app.config["MAX_CONTENT_LENGTH"] = 2 * 1024 * 1024  # File has to be less than 2MB
app.config.update(
    CELERY_CONFIG={
        "broker_url": "redis://localhost:6379",
        "result_backend": "redis://localhost:6379",
    }
)


def call_basher(sample_sheet):
    args = Namespace(
        input_spreadsheet=sample_sheet,
    )
    (
        mode_of_inheritance,
        sample_id,
        number_snps_imported,
        summary_snps_by_region,
        informative_snps_by_region,
        embryo_count_data_df,
        html_string,
    ) = excel_parser.main(args)
    return sample_id, html_string


class ChangeForm(FlaskForm):
    sample_sheet = FileField("Sample Sheet:")
    snp_array_files = MultipleFileField("SNP Array Files:")
    submit = SubmitField("Run BASHer")


@app.route("/", methods=["GET", "POST"])
def form(basher_state="initial"):
    SampleSheetUp = SampleSheetUpload()  # File Upload class - detailed below.
    SnpArrayUp = SnpArrayUpload()  # File Upload class - detailed below.
    chgForm = ChangeForm()
    chgDetail = dict()

    if request.method == "POST" and chgForm.validate_on_submit():
        sample_sheet = chgForm.sample_sheet.data
        x = SampleSheetUp.upload(sample_sheet)
        chgDetail["sample_sheet"] = x

        snp_array_files = chgForm.snp_array_files.data
        y = SnpArrayUp.upload(snp_array_files)
        chgDetail["snp_array_files"] = y
        basher_state = "started"

        sample_id, html_report = call_basher(x)
        timestr = datetime.now().strftime("%Y%m%d-%H%M%S")
        session["report_name"] = f"{sample_id}_{timestr}.html"
        session["report_path"] = os.path.join(
            app.config["UPLOAD_FOLDER"],
            session["report_name"],
        )
        with open(
            session["report_path"],
            "w",
        ) as f:
            f.write(html_report)

        return render_template(
            "index.html",
            form=chgForm,
            basher_state=basher_state,
            sample_sheet_name=sample_sheet.filename,
            snp_array_file_names=snp_array_files[0].filename,
            report_name=session["report_name"],
        )

    return render_template(
        "index.html",
        form=chgForm,
        basher_state=basher_state,
    )


class SampleSheetUpload:
    def allowed_files(self, file_name):
        allowed_extensions = [
            "xlsm",
            "xlsx",
        ]
        return (
            "." in file_name
            and file_name.rsplit(".", 1)[1].lower() in allowed_extensions
        )

    def upload(self, file):
        file_name = file.filename
        if file_name == "":
            return "NULL"

        elif file_name and self.allowed_files(file_name):
            secure_file_name = secure_filename(file_name)
            file.save(
                os.path.join(
                    app.config["UPLOAD_FOLDER"],
                    secure_file_name,
                )
            )

            return str(
                os.path.join(
                    app.config["UPLOAD_FOLDER"],
                    secure_file_name,
                )
            )


class SnpArrayUpload:
    def allowed_files(self, file_name):
        allowed_extensions = [
            "txt",
            "csv",
        ]
        return (
            "." in file_name
            and file_name.rsplit(".", 1)[1].lower() in allowed_extensions
        )

    def upload(self, file):
        # TODO add code to hand multiple files
        file_name = file[0].filename
        if file_name == "":
            return "NULL"

        elif file_name and self.allowed_files(file_name):
            secure_file_name = secure_filename(file_name)
            file[0].save(
                os.path.join(
                    app.config["UPLOAD_FOLDER"],
                    secure_file_name,
                )
            )

            return str(
                os.path.join(
                    app.config["UPLOAD_FOLDER"],
                    secure_file_name,
                )
            )


@app.route("/download")
def download():
    report_path = session["report_path"]
    file_name = session["report_name"]
    return send_file(
        report_path, download_name=file_name, mimetype="text/html", as_attachment=True
    )


if __name__ == "__main__":
    app.run(debug=True)
