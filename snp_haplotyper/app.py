from argparse import Namespace
from datetime import datetime
import excel_parser
from flask import Flask, render_template, Response, send_file, jsonify, request, session
from flask_wtf import FlaskForm
from flask_session import Session
import merge_array_files
import os
import pdfkit
import time
import random
from wtforms import FileField, SubmitField, MultipleFileField
from werkzeug.utils import secure_filename
import zipfile


import logging

log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
logger = logging.getLogger("BASHer_logger")

# To override the default severity of logging
logger.setLevel("DEBUG")


app = Flask(__name__)
app.config["UPLOAD_FOLDER"] = os.environ["UPLOAD_FOLDER"]
app.config["SECRET_KEY"] = "catchmeifyoucan"
app.config["SESSION_TYPE"] = "filesystem"
app.config["SESSION_FILE_DIR"] = os.environ["SESSION_FILE_DIR"]
app.config["UPLOAD_EXTENSIONS"] = [".txt", ".csv", ".xlsm", ".xlsx"]
app.config["MAX_CONTENT_LENGTH"] = 2 * 1024 * 1024  # File has to be less than 2MB
app.config.from_object(__name__)
Session(app)


def call_basher(sample_sheet, snp_array_file):
    args = Namespace(
        input_spreadsheet=sample_sheet,
        snp_array_file=snp_array_file,
    )
    (
        mode_of_inheritance,
        sample_id,
        number_snps_imported,
        summary_snps_by_region,
        informative_snps_by_region,
        embryo_count_data_df,
        html_string,
        pdf_string,
    ) = excel_parser.main(args)
    return sample_id, html_string, pdf_string


class ChangeForm(FlaskForm):
    sample_sheet = FileField("Sample Sheet:")
    snp_array_files = MultipleFileField("SNP Array Files:")
    submit = SubmitField("Run BASHer")


@app.route("/basher/", methods=["GET", "POST"])
def form(basher_state="initial"):
    SampleSheetUp = SampleSheetUpload()  # File Upload class - detailed below.
    SnpArrayUp = SnpArrayUpload()  # File Upload class - detailed below.
    chgForm = ChangeForm()
    chgDetail = dict()

    if request.method == "POST" and chgForm.validate_on_submit():
        # Use FileHandler() to log to a file
        session["timestr"] = datetime.now().strftime("%Y%m%d-%H%M%S")
        os.mkdir(os.path.join(app.config["UPLOAD_FOLDER"], session["timestr"]))
        file_handler = logging.FileHandler(f"/var/local/basher/logs/basher_error.log")
        formatter = logging.Formatter(log_format)
        file_handler.setFormatter(formatter)

        # Don't forget to add the file handler
        logger.addHandler(file_handler)

        sample_sheet = chgForm.sample_sheet.data
        input_sheet = SampleSheetUp.upload(sample_sheet)
        chgDetail["sample_sheet"] = input_sheet

        snp_array_files = chgForm.snp_array_files.data
        input_files = SnpArrayUp.upload(snp_array_files)
        # If multiple files are uploaded merge them into a single file, else just use the single file
        if len(input_files) > 1:
            # use removesuffix to remove .txt from the end of the file names in the list
            # and then join them together with a _ to make a new file name
            merged_file_name = (
                "_".join(
                    [os.path.basename(x).removesuffix(".txt") for x in input_files]
                )
                + ".txt"
            )
            df = merge_array_files.main(input_files)
            df.to_csv(
                os.path.join(
                    app.config["UPLOAD_FOLDER"], session["timestr"], merged_file_name
                ),
                sep="\t",
            )
            input_file = merged_file_name
        else:
            input_file = input_files[0]

        chgDetail["snp_array_files"] = input_file
        basher_state = "started"

        # Get the file names of the uploaded files
        input_sheet_basename = os.path.basename(input_sheet)
        input_file_basename = os.path.basename(input_file)

        # Create the paths to the uploaded files
        input_sheet_tmp_path = os.path.join(
            app.config["UPLOAD_FOLDER"], session["timestr"], input_sheet_basename
        )

        input_file_tmp_path = os.path.join(
            app.config["UPLOAD_FOLDER"], session["timestr"], input_file_basename
        )

        sample_id, html_report, pdf_report = call_basher(
            input_sheet_tmp_path, input_file_tmp_path
        )
        session["report_name"] = f'{sample_id}_{session["timestr"]}'
        session["report_path"] = os.path.join(
            app.config["UPLOAD_FOLDER"],
            session["report_name"],
        )
        with open(
            f'{session["report_path"]}.html',
            "w",
        ) as f:
            f.write(html_report)
        logger.info(
            f"Saved HTML report for {sample_id} at {session['report_path']}.html"
        )

        # Convert HTML report to PDF
        #        pdfkit.from_string(
        #            pdf_report,
        #            f'{session["report_path"]}.pdf',
        #        )
        #        logger.info(f"Saved PDF report for {sample_id}")

        return render_template(
            "index.html",
            form=chgForm,
            basher_state=basher_state,
            sample_sheet_name=sample_sheet.filename,
            snp_array_file_names=", ".join([x.filename for x in snp_array_files]),
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
                    session["timestr"],
                    secure_file_name,
                )
            )

            return str(
                os.path.join(
                    app.config["UPLOAD_FOLDER"],
                    session["timestr"],
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

    def upload(self, files):
        file_names = []

        for file in files:
            file_name = file.filename
            if file_name == "":
                return "NULL"

            elif file_name and self.allowed_files(file_name):
                secure_file_name = secure_filename(file_name)
                file.save(
                    os.path.join(
                        app.config["UPLOAD_FOLDER"],
                        session["timestr"],
                        secure_file_name,
                    )
                )
                file_names.append(
                    str(
                        os.path.join(
                            app.config["UPLOAD_FOLDER"],
                            session["timestr"],
                            secure_file_name,
                        )
                    )
                )

        return file_names


@app.route("/basher/download")
def download():
    html_file_name = f'{session["report_name"]}.html'
    #    pdf_file_name = f'{session["report_name"]}.pdf'

    # Create a zip file with the html and pdf reports
    with zipfile.ZipFile(f'{session["report_path"]}.zip', "w") as zipObj:
        zipObj.write(f'{session["report_path"]}.html', html_file_name)
    #       zipObj.write(f'{session["report_path"]}.pdf', pdf_file_name)

    logger.info(f'Saved zipped reports to {session["report_path"]}.zip')

    # Delete the html and pdf reports
    os.remove(f'{session["report_path"]}.html')
    #   os.remove(f'{session["report_path"]}.pdf')

    logger.info(f'Attempting to download {session["report_path"]}.html')

    return send_file(
        f'{session["report_path"]}.zip',
        download_name=f'{session["report_path"]}.zip',
        mimetype="text/html",
        as_attachment=True,
    )


if __name__ == "__main__":
    app.run(debug=True)
