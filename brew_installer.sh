#!/bin/sh

# This will ONLY work on macOS - High Sierra (10.13) or newer likely work, tested on macOS Monterey (12.5).

xcode-select --install
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

brew install python@3.10
brew install hdf5
pip3 install termcolor
