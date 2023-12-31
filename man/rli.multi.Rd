% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rli.R
\name{rli.multi}
\alias{rli.multi}
\title{Red List Index for multiple groups.}
\usage{
rli.multi(spData, tree = NULL, boot = FALSE, dd = FALSE, runs = 1000)
}
\arguments{
\item{spData}{A matrix with group names (first column) and species assessment categories for one or two points in time (remaining columns). Values can be text (EX, EW, RE, CR, EN, VU, NT, DD, LC) or numeric (0 for LC, 1 for NT, 2 for VU, 3 for EN, 4 for CR, 5 for RE/EW/EX).}

\item{tree}{A list of hclust or phylo objects, each corresponding to a tree per group (used when species are weighted by their unique contribution to phylogenetic or functional diversity).}

\item{boot}{If TRUE bootstrapping for statistical significance is performed on both values per date and the trend between dates.}

\item{dd}{bootstrap among all species (FALSE) or Data Deficient species only (TRUE).}

\item{runs}{Number of runs for bootstrapping}
}
\value{
A matrix with the RLI values and, if bootstrap is performed, their confidence limits and significance.
}
\description{
Calculates the Red List Index (RLI) for multiple groups of species.
}
\details{
The IUCN Red List Index (RLI) (Butchart et al. 2004, 2007) reflects overall changes in IUCN Red List status over time of a group of taxa.
The RLI uses weight scores based on the Red List status of each of the assessed species. These scores range from 0 (Least Concern) to 5 (Extinct/Extinct in the Wild).
Summing these scores across all species and relating them to the worst-case scenario, i.e. all species extinct, gives us an indication of how biodiversity is doing.
Each species weight can further be influenced by how much it uniquely contributes to the phylogenetic or functional diversity of the group (Cardoso et al. in prep.).
Importantly, the RLI is based on true improvements or deteriorations in the status of species, i.e. genuine changes. It excludes category changes resulting from, e.g., new knowledge (Butchart et al. 2007).
The RLI approach helps to develop a better understanding of which taxa, regions or ecosystems are declining or improving.
Juslen et al. (2016a, b) suggested the use of bootstrapping to search for statistical significance when comparing taxa or for trends in time of the index and this approach is here implemented.
}
\examples{
rliData <- matrix(c("LC","LC","EN","EN","EX","EX","LC","CR","CR","EX"), ncol = 2, byrow = TRUE)
colnames(rliData) <- c("2000", "2010")
rliData <- cbind(c("Arthropods","Arthropods","Birds","Birds","Birds"), rliData)
rli.multi(rliData[,1:2])
rli.multi(rliData[,1:2], boot = TRUE)
rli.multi(rliData)
rli.multi(rliData, boot = TRUE)
}
\references{
Butchart, S.H.M., Stattersfield, A.J., Bennun, L.A., Shutes, S.M., Akcakaya, H.R., Baillie, J.E.M., Stuart, S.N., Hilton-Taylor, C. & Mace, G.M. (2004) Measuring global trends in the status of biodiversity: Red List Indices for birds. PloS Biology, 2: 2294-2304.

Butchart, S.H.M., Akcakaya, H.R., Chanson, J., Baillie, J.E.M., Collen, B., Quader, S., Turner, W.R., Amin, R., Stuart, S.N. & Hilton-Taylor, C. (2007) Improvements to the Red List index. PloS One, 2: e140.

Juslen, A., Cardoso, P., Kullberg, J., Saari, S. & Kaila, L. (2016a) Trends of extinction risk for Lepidoptera in Finland: the first national Red List Index of butterflies and moths. Insect Conservation and Diversity, 9: 118-123.

Juslen, A., Pykala, J., Kuusela, S., Kaila, L., Kullberg, J., Mattila, J., Muona, J., Saari, S. & Cardoso, P. (2016b) Application of the Red List Index as an indicator of habitat change. Biodiversity and Conservation, 25: 569-585.
}
