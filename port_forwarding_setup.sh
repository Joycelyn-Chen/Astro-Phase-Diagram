# On local machine
ssh -L 2222:compute2:22 -L 8889:compute2:8889 portal
ssh -p 2222 localhost


# On remote server
jupyter notebook --no-browser --port=8889 --ip=0.0.0.0 

# Connect from local browser
http://localhost:8889


# Permanent solution 
nohup jupyter notebook --no-browser --port=8889 --ip=0.0.0.0 &
# (on remote server)