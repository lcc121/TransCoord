# TransCoord
TransCoord is a genome-guided transcriptome assembler for RNA-seq.

## Description
TransCoord is a sensitive genome-guided transcriptome assembler for RNA-seq data. It inclusively gathers all kinds of candidate transcripts into a two-phased linear programming model that aims both to minimize coverage deviation and reserve least transcripts under conserving the must. In this way, the outcome is a coordination of all candidates, instead of a union of all independently assembled parts.

## Installation
Because of the embedded Scallop usage, you have to install additional libraries used in [Scallop](https://github.com/Kingsford-Group/scallop/blob/master/README.md), including Boost and htslib.  
You also need to install a linear programming solver [Gurobi](https://www.gurobi.com/) ( free for academic license ).  
After install these dependencies, you can compile the source code of TransCoord.
### Install Boost and htslib
Same as in [Scallop](https://github.com/Kingsford-Group/scallop/blob/master/README.md).
### Install Gurobi
First download Gurobi optimizer from this [link](https://www.gurobi.com/downloads/gurobi-optimizer-eula/). Next follow the steps of this [link](https://www.gurobi.com/downloads/end-user-license-agreement-academic/) to apply a academic license for the free usage of Gurobi. Then click the *License ID* link to the *License Details* page and follow the instructions of Installation to setup the license file on your computer.
### Compile TransCoord

## Usage
### Quick Start

## Contack Information
If you have any questions or concerns, please feel free to contact <lichenchen121@outlook.com>.
