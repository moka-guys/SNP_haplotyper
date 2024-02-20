# pull official base image
FROM python:3.10.6-slim-buster

# set working directory
WORKDIR /usr/local/basher/snp_haplotyper

# Accept the commit hash as a build argument
ARG IMG_VERSIONED

# set environment variables
ENV IMG_VERSIONED ${IMG_VERSIONED}
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1
#ENV FLASK_DEBUG 1
#ENV FLASK_ENV development
ENV SCRIPT_NAME=""
ENV FLASK_APP /usr/local/basher/snp_haplotyper/app.py
ENV PYTHONPATH /usr/local/basher/snp_haplotyper
ENV UPLOAD_FOLDER /var/local/basher/uploads
ENV SESSION_FILE_DIR /var/local/basher/flask_sessions
# only for development

# add and install requirements
COPY ./requirements.txt .
RUN pip3 install -r requirements.txt
USER 0
RUN mkdir -p /var/local/basher/logs/
RUN mkdir -p /var/local/basher/uploads/
RUN mkdir -p /var/local/basher/flask_sessions/
RUN chmod 777 /var/local/basher/logs/
RUN chmod 777 /var/local/basher/uploads/
RUN chmod 777 /var/local/basher/flask_sessions/
USER $CONTAINER_USER_ID

RUN apt-get update && apt-get install -y wkhtmltopdf

# add app
COPY ["snp_haplotyper", "requirements.txt", "tests", "wsgi.py", \
    "pytest.ini", "docs", "gunicorn.conf.py", "./"]

COPY ["test_data/AffyID2rsid.txt", "../test_data/AffyID2rsid.txt"]

EXPOSE 5000

# run server
CMD ["gunicorn", "wsgi:app"]
