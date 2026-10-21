#!/usr/bin/env python3
"""List SRA runs for one or more BioProjects from the AWS Open Data copy of the
SRA metadata (s3://sra-pub-metadata-us-east-1), without NCBI or ENA access.

usage: sra_runs.py PRJNA1145481 [PRJEB9717 ...] > runs.tsv

Needs pyarrow. Each metadata file is read for the selected projects only, so
the whole 14 GB table is never downloaded.
"""
import os
import sys
from concurrent.futures import ThreadPoolExecutor

import pyarrow.fs as pafs
import pyarrow.parquet as pq

COLUMNS = ["acc", "bioproject", "library_name", "sample_name", "organism",
           "librarylayout", "libraryselection", "instrument", "mbases",
           "avgspotlen", "releasedate"]


def main(projects):
    proxy = os.environ.get("HTTPS_PROXY") or os.environ.get("https_proxy")
    fs = pafs.S3FileSystem(anonymous=True, region="us-east-1", proxy_options=proxy)
    selector = pafs.FileSelector("sra-pub-metadata-us-east-1/sra/metadata/")
    files = [f.path for f in fs.get_file_info(selector) if f.type == pafs.FileType.File]

    def read(path):
        return pq.read_table(path, columns=COLUMNS, filesystem=fs,
                             filters=[("bioproject", "in", projects)]).to_pylist()

    rows = []
    with ThreadPoolExecutor(4) as pool:
        for part in pool.map(read, files):
            rows.extend(part)
    print("\t".join(COLUMNS))
    for row in sorted(rows, key=lambda r: (r["bioproject"], r["library_name"] or "", r["acc"])):
        print("\t".join("" if row[c] is None else str(row[c]) for c in COLUMNS))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1:])
