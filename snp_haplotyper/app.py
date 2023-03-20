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

import logging

log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
logger = logging.getLogger("BASHer_logger")

# To override the default severity of logging
logger.setLevel("DEBUG")


app = Flask(__name__)
app.config["UPLOAD_FOLDER"] = os.environ["UPLOAD_FOLDER"]  # Set in docker-compose.yml
app.config["SECRET_KEY"] = "catchmeifyoucan"
app.config["UPLOAD_EXTENSIONS"] = [".txt", ".csv", ".xlsm", ".xlsx"]
app.config["MAX_CONTENT_LENGTH"] = 2 * 1024 * 1024  # File has to be less than 2MB


def call_basher(sample_sheet, snp_array_files):
    args = Namespace(
        input_spreadsheet=sample_sheet,
        input_files=snp_array_files,  # One or more SNP array files
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
        # Use FileHandler() to log to a file
        timestr = datetime.now().strftime("%Y%m%d-%H%M%S")
        file_handler = logging.FileHandler(f"{timestr}_basher_.log")
        formatter = logging.Formatter(log_format)
        file_handler.setFormatter(formatter)

        # Don't forget to add the file handler
        logger.addHandler(file_handler)

        sample_sheet = chgForm.sample_sheet.data
        input_sheet = SampleSheetUp.upload(sample_sheet)
        chgDetail["sample_sheet"] = input_sheet

        snp_array_files = chgForm.snp_array_files.data
        input_files = SnpArrayUp.upload(snp_array_files)
        chgDetail["snp_array_files"] = input_files
        basher_state = "started"

        # Get the file names of the uploaded files
        input_sheet_basename = os.path.basename(input_sheet)
        input_files_basenames = [os.path.basename(x) for x in input_files]

        # Create the paths to the uploaded files
        input_sheet_tmp_path = os.path.join(
            app.config["UPLOAD_FOLDER"], input_sheet_basename
        )

        input_files_tmp_paths = [
            os.path.join(app.config["UPLOAD_FOLDER"], x) for x in input_files_basenames
        ]

        sample_id, html_report = call_basher(
            input_sheet_tmp_path, input_files_tmp_paths
        )
        session["report_name"] = f"{sample_id}_{timestr}"
        session["report_path"] = os.path.join(
            app.config["UPLOAD_FOLDER"],
            session["report_name"],
        )
        with open(
            session["report_path"],
            "w",
        ) as f:
            f.write(html_report)
        logger.info(f"Saved HTML report for {sample_id}")

        return render_template(
            "index.html",
            form=chgForm,
            basher_state=basher_state,
            sample_sheet_name=sample_sheet.filename,
            snp_array_file_names=snp_array_files[0].filename,
            report_name=f'{session["report_name"]}.html',
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
    file_name = f'{session["report_name"]}.html'
    return send_file(
        report_path, download_name=file_name, mimetype="text/html", as_attachment=True
    )


if __name__ == "__main__":
    app.run(debug=True)
