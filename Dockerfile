FROM ubuntu
RUN mkdir /mwfitting
WORKDIR /mwfitting
ADD requirements.txt /mwfitting/requirements.txt
RUN apt-get update
RUN apt-get install -y python 
RUN apt-get install -y python-pip
 
# Install mwfitting dependencies
RUN pip install --upgrade pip
RUN /bin/bash -c "pip install -r /mwfitting/requirements.txt"

RUN useradd -d /home/docker -m docker
USER docker
CMD [ "/bin/bash" ]
