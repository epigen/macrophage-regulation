#!/bin/bash

env_name="$1"

if [ -z "$env_name" ]; then
    echo "Usage: $0 <env_name>"
    exit 1
fi

conda env export --no-build --from-history -n "$env_name" > "./$env_name.fromHistory.yaml"
conda env export --no-build -n "$env_name" > "./$env_name.noBuild.yaml"
conda env export -n "$env_name" > "./$env_name.all.yaml"

# Removing the last line from each file containing the local path of the environment
sed -i '$ d' "./$env_name.fromHistory.yaml"
sed -i '$ d' "./$env_name.noBuild.yaml"
sed -i '$ d' "./$env_name.all.yaml"
