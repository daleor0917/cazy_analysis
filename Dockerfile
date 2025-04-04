FROM gitpod/workspace-full

# Instalar R
RUN sudo apt update && sudo apt install -y r-base

# Crear un directorio en HOME con permisos para librerías de R
RUN mkdir -p /home/gitpod/Rlibs

# Instalar paquetes de R en ese directorio
RUN R -e 'install.packages(c("shiny", "reticulate"), repos="https://cloud.r-project.org/", lib="/home/gitpod/Rlibs")'

# Establecer la variable de entorno para que R siempre use este directorio
ENV R_LIBS_USER="/home/gitpod/Rlibs"

