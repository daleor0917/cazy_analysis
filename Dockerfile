# Usar Python base
FROM python:3.10-slim

# Instalar dependencias mínimas del sistema
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    mafft \
    iqtree \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# Descargar y mover el binario precompilado de ModelTest-NG
RUN wget https://github.com/ddarriba/modeltest/releases/download/v0.1.7/modeltest-ng-0.1.7-linux.tar.gz && \
    tar -xzf modeltest-ng-0.1.7-linux.tar.gz && \
    mv modeltest-ng /usr/local/bin/ && \
    chmod +x /usr/local/bin/modeltest-ng && \
    rm -rf modeltest-ng-0.1.7-linux.tar.gz

# Instalar paquetes de Python necesarios
RUN pip install --no-cache-dir ete3 PyQt5

# Crear directorio de trabajo
WORKDIR /app

# Copiar archivos del proyecto
COPY . /app