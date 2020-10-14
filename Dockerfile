FROM python:3.7.9-buster


RUN apt-key adv --keyserver keys.gnupg.net --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF' && \
    echo 'deb http://mirrors.aliyun.com/debian/ buster main non-free contrib' > /etc/apt/sources.list && \
    echo 'deb-src http://mirrors.aliyun.com/debian/ buster main non-free contrib' >> /etc/apt/sources.list && \
    echo 'deb http://mirrors.aliyun.com/debian-security buster/updates main' >> /etc/apt/sources.list && \
    echo 'deb-src http://mirrors.aliyun.com/debian-security buster/updates main' >> /etc/apt/sources.list && \
    echo 'deb http://mirrors.aliyun.com/debian/ buster-updates main non-free contrib' >> /etc/apt/sources.list && \
    echo 'deb-src http://mirrors.aliyun.com/debian/ buster-updates main non-free contrib' >> /etc/apt/sources.list && \
    echo 'deb http://mirrors.aliyun.com/debian/ buster-backports main non-free contrib' >> /etc/apt/sources.list && \
    echo 'deb-src http://mirrors.aliyun.com/debian/ buster-backports main non-free contrib' >> /etc/apt/sources.list && \
    echo 'deb http://cloud.r-project.org/bin/linux/debian buster-cran40/' >> /etc/apt/sources.list && \
    apt update && apt install -y samtools bedtools r-base


RUN pip install -i https://mirrors.aliyun.com/pypi/simple click peakutils pysam pybedtools matplotlib loguru sklearn

RUN echo 'options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))' > ~/.Rprofile && \
    Rscript -e 'install.packages(c("ggplot2", "reshape2", "dplyr"))'

COPY ./ /opt

ENTRYPOINT [ "python", "/opt/main.py" ]