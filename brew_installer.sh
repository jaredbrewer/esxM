#!/bin/sh

# This will ONLY work on Macs - macOS Mojave (10.14) or High Sierra (10.13) preferred, untested on other versions.

xcode-select --install
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

brew tap homebrew/science
brew install kallisto
brew install python@3.10
brew install --cask R
brew install --cask rstudio
brew install rpy2
pip3 install BioPython
pip3 install termcolor
