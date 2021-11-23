FROM us.gcr.io/anvil-gcr-public/anvil-rstudio-bioconductor-devel:3.14.2

RUN apt update; apt install -y nvidia-cuda-toolkit
RUN wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-ubuntu2004.pin
RUN sudo mv cuda-ubuntu2004.pin /etc/apt/preferences.d/cuda-repository-pin-600
RUN wget https://developer.download.nvidia.com/compute/cuda/11.4.2/local_installers/cuda-repo-ubuntu2004-11-4-local_11.4.2-470.57.02-1_amd64.deb
RUN sudo dpkg -i cuda-repo-ubuntu2004-11-4-local_11.4.2-470.57.02-1_amd64.deb
RUN rm cuda-repo-ubuntu2004-11-4-local_11.4.2-470.57.02-1_amd64.deb
RUN sudo apt-key add /var/cuda-repo-ubuntu2004-11-4-local/7fa2af80.pub
RUN sudo apt-get update
RUN sudo apt-get -y install cuda

COPY cudnn-11.4-linux-x64-v8.2.4.15.tgz /tmp
RUN tar zxf /tmp/cudnn-11.4-linux-x64-v8.2.4.15.tgz
RUN cp cuda/include/cudnn.h /usr/lib/cuda/include/
RUN cp cuda/lib64/libcudnn* /usr/lib/cuda/lib64/
RUN chmod a+r /usr/lib/cuda/include/cudnn.h /usr/lib/cuda/lib64/libcudnn*

ENV LD_LIBRARY_PATH="/usr/lib/cuda/lib64:$LD_LIBRARY_PATH"

RUN pip3 install tensorflow==2.2.0
