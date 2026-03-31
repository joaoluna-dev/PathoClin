FROM python:3.10-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential python3-dev perl bcftools samtools tabix \
    libpango-1.0-0 libpangoft2-1.0-0 libjpeg-dev libopenjp2-7-dev \
    libffi-dev shared-mime-info libcairo2 zlib1g-dev libbz2-dev \
    liblzma-dev libcurl4-openssl-dev libssl-dev fonts-liberation \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements.txt /app/
RUN pip install --no-cache-dir -r requirements.txt

COPY . /app/

RUN chmod +x /app/run.sh
RUN chmod +x /app/clean_data.sh

EXPOSE 8501

CMD ["./run.sh"]