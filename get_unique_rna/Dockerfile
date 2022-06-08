FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:9a7d-main

# Its easy to build binaries from source that you can later reference as
# subprocesses within your workflow.

COPY ./Alloc_RNA /root/bin/Alloc_RNA
RUN chmod +x -R /root/bin/Alloc_RNA

# Make sure the c++ library is available
RUN export LD_LIBRARY_PATH=/root/bin/Alloc_RNA:$LD_LIBRARY_PATH && ldconfig

RUN pip install biopython

ENV PYTHONPATH "${PYTHONPATH}:/root/wf"

# STOP HERE:
# The following lines are needed to ensure your build environement works
# correctly with latch.
COPY wf /root/wf
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
RUN python3 -m pip install --upgrade latch
WORKDIR /root
