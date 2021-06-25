FROM continuumio/anaconda3

WORKDIR /opt/scvr

ENV SHELL bash

ADD ./requirements.txt /tmp/requirements.txt
# Add our code
ADD . /opt/scvr

# Install dependencies
RUN pip install -r /tmp/requirements.txt

#CMD conda run -n env_stream 
CMD gunicorn --pythonpath dash_app/apps/singlecell-vr-api app:server --timeout 300
