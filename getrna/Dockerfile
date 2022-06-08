FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/wf-base:fbe8-main

# You can use local data to construct your workflow image.  Here we copy a
# pre-indexed reference to a path that our workflow can reference.
COPY data /root/data

# Include the Alloc_RNA executable in the image.
COPY /bin/Alloc_RNA /root/bin

# STOP HERE:
# The following lines are needed to ensure your build environement works
# correctly with latch.
COPY wf /root/wf
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
RUN python3 -m pip install --upgrade latch
WORKDIR /root
