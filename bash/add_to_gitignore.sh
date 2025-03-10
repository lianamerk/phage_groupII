#!/bin/bash

# Define the size threshold in bytes (100MB in this example)
size_threshold=$((100 * 1024 * 1024)) # 100 MB

# Find files exceeding the size threshold and append them to .gitignore
find . -type f -exec stat -c "%s %n" {} + | awk -v size_threshold=$size_threshold '$1 > size_threshold { print $2 }' >> .gitignore
