#!/bin/bash

# Instalar SDKMAN si no está instalado
if [ ! -d "$HOME/.sdkman" ]; then
  curl -s "https://get.sdkman.io" | bash
  source "$HOME/.sdkman/bin/sdkman-init.sh"
else
  source "$HOME/.sdkman/bin/sdkman-init.sh"
fi

# Instalar Java 17 usando SDKMAN
sdk install java 17.0.8-tem
sdk default java 17.0.8-tem

# Verificar la instalación de Java
java -version

# Instalar Nextflow
curl -s https://get.nextflow.io | bash
mkdir -p ~/bin
mv nextflow ~/bin/
chmod +x ~/bin/nextflow

# Asegurarse de que ~/bin esté en el PATH
gp env PATH="$HOME/bin:$PATH"