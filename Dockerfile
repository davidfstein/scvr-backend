FROM continuumio/anaconda3

ADD ./requirements.txt /tmp/requirements.txt

# Add our code
ADD . /opt/scvr
WORKDIR /opt/scvr

ENV SHELL bash

RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge 
  
#RUN apt-get update && apt-get install gsl-bin libgsl0-dev -y && apt-get clean
#RUN conda create -n env_stream python stream=1.0
RUN conda create -n env_stream python stream=1.0

SHELL ["conda", "run", "-n", "env_stream", "/bin/bash", "-c"]

# Install dependencies
RUN pip install -r /tmp/requirements.txt

CMD conda run -n env_stream gunicorn --pythonpath dash_app/apps/singlecell-vr-api app:server --timeout 300
