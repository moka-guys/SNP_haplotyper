# file gunicorn.conf.py
# coding=utf-8

import os
import multiprocessing

_LOGS = '/var/local/basher/logs'
errorlog = os.path.join(_LOGS, 'basher_error.log')
accesslog = os.path.join(_LOGS, 'basher_access.log')
loglevel = 'info'
bind = '0.0.0.0:5000'
workers = 2 # multiprocessing.cpu_count() * 2 + 1
timeout = 30  # timeout 30 seconds
keepalive = 60 * 60  # keep connections alive for 1 hour
capture_output = True
