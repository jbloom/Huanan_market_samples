#!/bin/bash

# generate track changes LaTex

latexdiff \
    --exclude-safecmd=citep \
    --exclude-safecmd=cite \
    paper_original_submission.tex \
    paper.tex \
    > paper_track_changes.tex
