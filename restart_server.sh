#!/bin/bash

# URL to check
URL="http://10.189.213.201/basher/"

# Command to restart the server (replace this with the actual command)
# TODO check the file path for this command on the server
RESTART_COMMAND="ansible-playbook -i graeme/deployment/inventories/production graeme/deployment/playbooks/basher.yml"

# Perform the HTTP request using curl & check for successful exit code
if ! curl --output /dev/null --silent --head --fail "$URL"; then
  logger "The BASHer server is not running. Attempting to restart..."
  $RESTART_COMMAND
  logger "Restart command sent."
else
  logger "The BASHer server is running."
fi
