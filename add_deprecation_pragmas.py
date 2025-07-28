#!/usr/bin/env python3
import os

# Configuration
old_version_dir = "include/legacy"  # Path to old headers
project_name = "DROMON"  # Project name for include path

# Ensure directory exists
if not os.path.exists(old_version_dir):
    print(f"Directory {old_version_dir} does not exist.")
    exit(1)

# Process headers
for filename in os.listdir(old_version_dir):
    if filename.endswith(".h") or filename.endswith(".hpp"):
        file_path = os.path.join(old_version_dir, filename)

        # Read original content
        with open(file_path, "r") as f:
            original_content = f.readlines()

        # Prepare new content
        pragma_lines = [
            "#pragma once\n",
            f'#pragma message("Warning: This header ({filename}) is deprecated. '
            f'Use <{project_name}/{filename}> instead.")\n\n',
            f'#include "{project_name}/{filename}"\n\n',
            "// ----------------- Original content below (commented out) -----------------\n"
        ]
        commented_content = [f"// {line}" for line in original_content]

        # Combine and write back
        with open(file_path, "w") as f:
            f.writelines(pragma_lines + commented_content)

        print(f"Updated: {file_path}")
