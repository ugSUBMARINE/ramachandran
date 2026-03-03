import argparse
import os
import time
from datetime import datetime


def cleanup_uploads(uploads_dir="uploads", max_age_hours=24, logger=None):
    """
    Deletes files in the specified directory that are older than max_age_hours.
    """

    def log_info(message):
        if logger:
            logger.info(message)
        else:
            print(message)

    def log_error(message):
        if logger:
            logger.error(message)
        else:
            print(message)

    if not os.path.exists(uploads_dir):
        log_info(f"Directory {uploads_dir} does not exist. Skipping cleanup.")
        return {"deleted": 0, "errors": 0}

    now = time.time()
    max_age_seconds = max_age_hours * 3600

    count = 0
    errors = 0

    log_info(f"Starting cleanup of {uploads_dir} (older than {max_age_hours} hours)...")

    for filename in os.listdir(uploads_dir):
        # Skip hidden files like .DS_Store or .gitkeep
        if filename.startswith("."):
            continue

        filepath = os.path.join(uploads_dir, filename)

        # Only process files
        if not os.path.isfile(filepath):
            continue

        file_age = os.path.getmtime(filepath)

        if (now - file_age) > max_age_seconds:
            try:
                os.remove(filepath)
                log_info(f"Deleted: {filename}")
                count += 1
            except Exception as e:
                log_error(f"Error deleting {filename}: {e}")
                errors += 1

    log_info(f"Cleanup finished. Deleted {count} files. Errors: {errors}.")
    return {"deleted": count, "errors": errors}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cleanup old files in the uploads directory.")
    parser.add_argument("--dir", default="uploads", help="Directory to clean up (default: uploads)")
    parser.add_argument(
        "--age",
        type=float,
        default=24,
        help="Maximum age of files in hours (default: 24)",
    )

    args = parser.parse_args()
    cleanup_uploads(args.dir, args.age)
