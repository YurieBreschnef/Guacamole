#!/bin/bash
echo "bash: uncompressing output data.."
unzip data.zip 
echo "bash: data.zip was created"
yes | rm data.zip
