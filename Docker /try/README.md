
#
chmod +x Analysis.sh
build
mkdir output
docker run --rm \
  -v /Users/cecilia.carmignoto/Documents/GitHub//Micro-C_pipeline/try/output:/output \
  --platform linux/amd64 rundirectly .