"""
Find a new RNA sequence not present in a given set of sequences. Allow users to upload a transcriptome or choose from a set of pre-loaded transcriptomes.
"""

import subprocess
from pathlib import Path

from latch import small_task, workflow
from latch.types import LatchFile
from enum import Enum
from typing import List, Dict, Tuple, Union, Optional

import import_fasta

class transcriptome(Enum):
    """
    Enum for transcriptome options
    """
    human = "human"
    mouse = "mouse"
    none = "none"

def get_transcriptome_file(transcriptome: transcriptome) -> Path:
    """
    Get the transcriptome file for the given transcriptome option
    """
    if transcriptome == transcriptome.human:
        return Path("/root/data/HumanTranscriptome.txt.dat")
    elif transcriptome == transcriptome.mouse:
        return Path("/root/data/MouseTranscriptome.txt.dat")

def GetTranscriptomeFromFile(transcriptome_file: Path) -> Path:
    #if the transcriptome is a text or dat file, return it's path
    if transcriptome_file.suffix == ".txt" or transcriptome_file.suffix == ".dat":
        return transcriptome_file
    #if the transcriptome is a FASTA file, convert it to a text file and return the path
    import_fasta.load_from_file(transcriptome_file, transcriptome_file.with_suffix(".txt"))
    
    return transcriptome_file.with_suffix(".txt")
    

@small_task
def get_RNA(length: int, transcriptome: transcriptome, transcriptome_file: Optional[LatchFile], additional_sequences: "list[str]" = []) -> str:

    transcriptome_file_path = ""

    # Are we using a pre-loaded transcriptome?
    if transcriptome == transcriptome.none:
        # If not, we need to upload a transcriptome file
        if transcriptome_file is None:
            raise ValueError("No transcriptome file provided")
        else:
            transcriptome_file_path = GetTranscriptomeFromFile(Path(transcriptome_file).resolve())
    else:
        # If we are using a pre-loaded transcriptome, we need to get the file path
        transcriptome_file_path = get_transcriptome_file(transcriptome)

    Alloc_cmd = [
        "/root/bin/Alloc_RNA",
        str(length),
        transcriptome_file_path,
        str(len(additional_sequences)),
       " ".join("\"{}\"".format(s) for s in additional_sequences)
    ]

    output = subprocess.run(Alloc_cmd, capture_output=True)

    return output.stdout.decode("utf-8")


@workflow
def get_unique_RNA(length: int, transcriptome: transcriptome, transcriptome_file: Optional[LatchFile], additional_sequences: "list[str]" = []) -> str:
    """Description...

    Get unique RNA
    --------------
    Get unique RNA from a given transcriptome. User may upload additional files or choose from a pre-loaded set of transcriptomes.

    __metadata__:
        display_name: Get unique RNA from a given transcriptome
        author:
            name: Cameron Kroll
            email: cameron.kroll@gmail.com
            github: github.com/camelcasecam
        repository: github.com/camelcasecam/alloc_rna
        license:
            id: MIT

    Args:

        length:
          Length of the RNA to be generated (must be greater than 10).

          __metadata__:
            display_name: Length of the RNA to be generated
        
        transcriptome:
            Transcriptome to check against when generating the RNA.

            __metadata__:
                display_name: Transcriptome to check against
        
        transcriptome_file:
            File containing the transcriptome to check against when generating the RNA (optional).

            __metadata__:
                display_name: Transcriptome file to check against (optional)
        
        additional_sequences:
            Additional sequences to check against when generating the RNA (optional).

            __metadata__:
                display_name: Additional sequences to check against (optional)
    """
    return get_RNA(length, transcriptome, transcriptome_file, additional_sequences)
