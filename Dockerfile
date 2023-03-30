# pull official base image
FROM python:3.10.6-slim-buster

# set working directory
WORKDIR /usr/local/basher

# set environment variables
#ENV PYTHONDONTWRITEBYTECODE 1
#ENV PYTHONUNBUFFERED 1
#ENV FLASK_DEBUG 1
#ENV FLASK_ENV development
ENV FLASK_APP /usr/local/basher/snp_haplotyper/app.py
ENV PYTHONPATH /usr/local/basher/snp_haplotyper
ENV UPLOAD_FOLDER /var/local/basher/uploads
# only for development

# add and install requirements
COPY ./requirements.txt .
RUN pip3 install -r requirements.txt
#RUN pip install debugpy

# add app
COPY . .

EXPOSE 5000

# run server
CMD ["gunicorn", "snp_haplotyper:app:app"]
