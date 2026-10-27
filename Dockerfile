# Sandbox for Inspect runs of the silicyte_bio_evals task: the same Python and packages the CI graders use. No network at run time.
FROM python:3.12-slim
RUN python -m pip install --no-cache-dir numpy==2.1.3 scipy==1.14.1 gemmi==0.6.5
WORKDIR /root
