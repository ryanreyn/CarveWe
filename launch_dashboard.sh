#!/bin/bash
# Launch CarveWe Dashboard from here
cd "$(dirname "$0")/carvewe-dashboard"
Rscript -e "shiny::runApp(port=3838, launch.browser=TRUE)"