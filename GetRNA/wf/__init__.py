"""
Find a new RNA sequence not present in a given set of sequences. Allow users to upload a transcriptome or choose from a set of pre-loaded transcriptomes.
"""

import subprocess
from pathlib import Path

from latch import small_task, workflow
from latch.types import LatchFile


@small_task
def assembly_task(length: int, transcriptome: LatchFile, transcriptome_option: str) -> LatchFile:

    # A reference to our output.
    sam_file = Path("covid_assembly.sam").resolve()

    Alloc_cmd = [
        "Alloc_RNA",
        str(length),
    ]

    subprocess.run(Alloc_cmd)

    return LatchFile(str(sam_file), "latch:///covid_assembly.sam")


@small_task
def sort_bam_task(sam: LatchFile) -> LatchFile:

    bam_file = Path("covid_sorted.bam").resolve()

    _samtools_sort_cmd = [
        "samtools",
        "sort",
        "-o",
        str(bam_file),
        "-O",
        "bam",
        sam.local_path,
    ]

    subprocess.run(_samtools_sort_cmd)

    return LatchFile(str(bam_file), "latch:///covid_sorted.bam")


@workflow
def get_unique_RNA(length: int, ) -> LatchFile:
    """Description...

    Get unique RNA
    --------------
    Get unique RNA from a given transcriptome. User may upload additional files or choose from a pre-loaded set of transcriptomes.

    __metadata__:
        display_name: Get unique RNA from a given transcriptome
        author:
            name: Cameron Kroll
            email: cameron.kroll@gmail.com
            github: github.com/CamelCaseCam
        repository: 
        license:
            id: MIT

    Args:

        read1:
          Paired-end read 1 file to be assembled.

          __metadata__:
            display_name: Read1

        read2:
          Paired-end read 2 file to be assembled.

          __metadata__:
            display_name: Read2
    """
    sam = assembly_task(read1=read1, read2=read2)
    return sort_bam_task(sam=sam)
