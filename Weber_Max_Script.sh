# Create empty directory for raw data
mkdir 01_raw

# Copy raw data (3 files) from server (server address may need to be adjusted for replication by a different user)
scp inf-27-2025@bioinf-serv2.cob.lu.se:/home2/resources/binp28/Data/ProjTaxa* ./01_raw/

# Unpack Data
gunzip ./01_raw/*.gz
