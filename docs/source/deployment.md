# Deployment

## Get BASHer from GitHub Repo

BASHer can be deployed on web server using docker image and ansible playbook. Create and push docker image using github repo https://github.com/moka-guys/SNP_haplotyper. Update the docker image name in `basher.yml` and deployment can be done on the web server using the command line 
```
ansible-playbook -i inventories/production playbooks/basher.yml
```
