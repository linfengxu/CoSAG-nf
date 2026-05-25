FROM python:3.8.6-slim
LABEL maintainer="Linfeng Xu"
LABEL description="Bio Python 3.8 environment"

# Debian Buster 已 EOL，替换为归档源
RUN sed -i 's|http://deb.debian.org/debian|http://archive.debian.org/debian|g' /etc/apt/sources.list \
    && sed -i 's|http://security.debian.org/debian-security|http://archive.debian.org/debian-security|g' /etc/apt/sources.list \
    && sed -i '/buster-updates/d' /etc/apt/sources.list

RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    libfreetype6-dev \
    libpng-dev \
    pkg-config \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --upgrade pip==21.0.1

RUN pip install --no-cache-dir \
    -i https://pypi.tuna.tsinghua.edu.cn/simple \
    beautifulsoup4==4.12.3 \
    biopython==1.83 \
    contourpy==1.1.1 \
    cycler==0.12.1 \
    fonttools==4.57.0 \
    importlib-resources==6.4.5 \
    joblib==1.4.2 \
    kiwisolver==1.4.7 \
    matplotlib==3.7.5 \
    numpy==1.24.4 \
    packaging==25.0 \
    pandas==2.0.3 \
    pillow==10.4.0 \
    pyparsing==3.1.4 \
    python-dateutil==2.9.0.post0 \
    pytz==2025.2 \
    PyYAML==6.0.2 \
    scikit-learn==1.3.2 \
    scipy==1.10.1 \
    seaborn==0.13.2 \
    six==1.17.0 \
    soupsieve==2.6 \
    threadpoolctl==3.5.0 \
    tzdata==2025.2 \
    zipp==3.20.2

WORKDIR /workspace
CMD ["python3"]
