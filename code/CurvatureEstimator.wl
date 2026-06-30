(* ::Package:: *)

(*
  PathAnisotropyCurvature / CurvatureEstimator.wl

  Updated clean core code for the shortest-path anisotropy diagnostic.

  Main estimator:
    C_log(p, r_g) = cubic mean deviation of Log[number of shortest paths]
    from p to vertices on the graph-distance shell r_g.

  Important convention:
    The main black-hole embedding benchmarks use the calibrated epsilon
    geometric graph construction:
        eps = distanceMatrixRadius[pts, k]
        g   = makeGeometricGraph[pts, eps]

    The diagnostic is a graph observable, not a graph-independent continuum
    scalar. Keep the graph construction and matched-control convention fixed
    when comparing results.
*)

BeginPackage["PathAnisotropyCurvature`"]; 

(* Numerical helpers *)
validRealNumberQ::usage = "validRealNumberQ[x] returns True when x is a finite real numeric value.";
safeStd::usage = "safeStd[x] returns StandardDeviation[x], or 0 for fewer than two values.";
safeCorrelation::usage = "safeCorrelation[x, y] computes Pearson correlation after removing invalid numeric pairs.";
rankData::usage = "rankData[x] returns average ranks, with ties assigned their mean rank.";
spearmanCorr::usage = "spearmanCorr[x, y] computes a Spearman rank correlation with basic validity checks.";
cubicMeanDeviation::usage = "cubicMeanDeviation[x] computes (Mean[Abs[x-Mean[x]]^3])^(1/3).";

(* Graph construction *)
distanceMatrixRadius::usage = "distanceMatrixRadius[pts, k, factor] returns factor times the median distance to the k-th nearest neighbor.";
makeGeometricGraph::usage = "makeGeometricGraph[pts, eps] builds an unweighted epsilon geometric graph from point coordinates.";
makeKNNGraph::usage = "makeKNNGraph[points, k] builds an unweighted k-nearest-neighbor graph from point coordinates.";

(* Point clouds and datasets *)
flatDiskPoints::usage = "flatDiskPoints[n] samples n points in a flat annulus/disk.";
hyperboloidPoints::usage = "hyperboloidPoints[n] samples n points on a simple hyperboloid-like benchmark surface.";
flammPoints::usage = "flammPoints[n, M] samples n points on a Flamm paraboloid z = 2 Sqrt[2 M (r - 2 M)].";
buildFlatDataset::usage = "buildFlatDataset[n, k] builds a flat epsilon-graph dataset association.";
buildHyperboloidDataset::usage = "buildHyperboloidDataset[n, k] builds a hyperboloid epsilon-graph benchmark dataset association.";
buildFlammDataset::usage = "buildFlammDataset[n, M, k] builds a Flamm/Schwarzschild epsilon-graph benchmark dataset association.";
buildMatchedFlatFromFlamm3::usage = "buildMatchedFlatFromFlamm3[flammData, k] makes a matched-flat kNN control by removing the Flamm height coordinate.";
buildMatchedFlatGeometricFromFlamm::usage = "buildMatchedFlatGeometricFromFlamm[flammData, k] makes a matched-flat epsilon-graph control by removing the Flamm height coordinate.";

(* Black-hole metric helpers *)
schwarzschildK::usage = "schwarzschildK[r, M] gives 48 M^2/r^6.";
rnF::usage = "rnF[M, Q][r] is the Reissner-Nordstrom lapse function.";
rnKretschmann::usage = "rnKretschmann[M, Q][r] is the Reissner-Nordstrom Kretschmann scalar.";
bardeenF::usage = "bardeenF[M, g][r] is the Bardeen lapse function.";
haywardF::usage = "haywardF[M, ell][r] is the Hayward lapse function.";
metricKretschmannFromF::usage = "metricKretschmannFromF[fexpr, r] computes K = f''^2 + (2 f'/r)^2 + 4(1-f)^2/r^4 for a symbolic expression fexpr.";
metricKretschmannFunctionFromF::usage = "metricKretschmannFunctionFromF[f] returns a numerical K(r) function from a numerical f[r].";
embeddingZValuesFromF::usage = "embeddingZValuesFromF[f, rs] numerically integrates z'(r)=Sqrt[1/f(r)-1] over the supplied radii.";
buildStaticSphericalFromExistingSampling::usage = "buildStaticSphericalFromExistingSampling[baseData, f, targetK, k, label] builds an epsilon-graph embedded surface using the r,phi samples from baseData.";
buildStaticSphericalFromFlammSampling::usage = "buildStaticSphericalFromFlammSampling[n, samplingMass, f, targetK, k, label] builds a static spherical embedding from Flamm-style sampling.";
buildRNDataset::usage = "buildRNDataset[n, M, Q, k] builds a Reissner-Nordstrom embedding dataset.";
buildBardeenDataset::usage = "buildBardeenDataset[n, M, g, k] builds a Bardeen embedding dataset.";
buildHaywardDataset::usage = "buildHaywardDataset[n, M, ell, k] builds a Hayward embedding dataset.";
buildBlackHoleFigureDatasets::usage = "buildBlackHoleFigureDatasets[n, M, k, Q, g, ell] returns an association of Flamm, RN, Bardeen, and Hayward datasets aligned to the same Flamm sampling. Zero limits are exact Flamm.";

(* Plotting utilities *)
plotDataset2D::usage = "plotDataset2D[data, label] plots the first two coordinates of a dataset.";
plotDataset3D::usage = "plotDataset3D[data, label] plots the first three coordinates of a dataset.";
shellAndBallVertices::usage = "shellAndBallVertices[data, center, radius] returns center, shell, and ball vertex sets.";
plotGeometryShell2D::usage = "plotGeometryShell2D[data, center, radius, label] colors center, ball, and shell vertices in 2D.";
plotGeometryShell3D::usage = "plotGeometryShell3D[data, center, radius, label] colors center, ball, and shell vertices in 3D.";
coloredGeometryPlot::usage = "coloredGeometryPlot[data, center, radius, label] chooses a 2D or 3D colored shell/ball plot automatically.";
estimatorValuesByVertex::usage = "estimatorValuesByVertex[data, res, estimator] returns a list indexed by vertex, containing estimator values or Missing.";
allGoodEstimatorValues::usage = "allGoodEstimatorValues[{{data,res},...}] returns all finite estimator values across panels.";
redBlueCF::usage = "redBlueCF[t] is a red-to-blue color function for t in [0,1].";
blueRedCF::usage = "blueRedCF[t] is a blue-to-red color function for t in [0,1].";
estimatorColoredGraph3DSharedScale::usage = "estimatorColoredGraph3DSharedScale[data, res, range, title] plots a 3D graph colored by LogCMD using a shared color scale.";
makeBlackHoleFamilyPanel::usage = "makeBlackHoleFamilyPanel[pairs, labels, range, colorFunction] returns a 2x2 shared-scale colored black-hole family panel.";

(* Estimator *)
precomputeAdjacency::usage = "precomputeAdjacency[g] returns an Association vertex -> adjacency list.";
shortestPathCountsFromAdj::usage = "shortestPathCountsFromAdj[adj, source] returns graph distances and shortest-path counts from source.";
pathAnisotropyEstimator::usage = "pathAnisotropyEstimator[g, p, radius] evaluates shortest-path anisotropy around vertex p.";
evaluateDataset::usage = "evaluateDataset[data, radius] evaluates the estimator for every vertex in a dataset.";
rowAssoc::usage = "rowAssoc[row] converts Dataset rows/rules/associations into an Association where possible.";
getRowValue::usage = "getRowValue[row, key] safely extracts a row value from string or symbol keys.";
cleanRows::usage = "cleanRows[res, estimator] converts Dataset/list rows into clean associations with valid estimator and TargetK values.";
estimatorStats::usage = "estimatorStats[res, estimator] computes Pearson and Spearman statistics against TargetK.";

(* Scans and binning *)
radiusScan::usage = "radiusScan[data, radii] evaluates estimatorStats for each graph radius.";
fastRadiusScan::usage = "fastRadiusScan[data, radii] is an alias of radiusScan.";
multiSeedRadiusScan::usage = "multiSeedRadiusScan[n, seeds, radii, k] runs Flamm benchmark scans across random seeds.";
summarizeSeedScan::usage = "summarizeSeedScan[scan] summarizes multi-seed scan results by graph radius.";
assignRadialBins::usage = "assignRadialBins[rows, nbins] adds a Bin key to rows based on their radial coordinate r.";
binnedFlammComparison::usage = "binnedFlammComparison[data, res, estimator, nbins] bins Flamm results radially and averages K and estimator values.";
radialBinComparisonFast2::usage = "radialBinComparisonFast2[data, res, estimator, nbins] bins generic results radially using cleaned rows.";
radialBinComparisonNoClean::usage = "radialBinComparisonNoClean[data, res, estimator, nbins] bins generic results radially without requiring TargetK.";
corrFromBinned::usage = "corrFromBinned[binned] computes Corr(rMean, EstimatorMean).";
curvatureScoreBins::usage = "curvatureScoreBins[binned] adds CurvatureScoreMean = -EstimatorMean, i.e. the reporting convention C_curv = -C_log.";
canonicalCurvatureBinnedCorrelations::usage = "canonicalCurvatureBinnedCorrelations[binned] computes Corr(rMean, C_curv) and Corr(Log[KMean], C_curv), with C_curv = -C_log.";
commonRadialBinnedComparison::usage = "commonRadialBinnedComparison[data1,res1,data2,res2,estimator,nbins] bins two datasets on the same overlapping radial domain.";
plotCommonRadialProfiles::usage = "plotCommonRadialProfiles[commonBinned,label1,label2,title,file] plots common-bin radial profiles and optionally exports them.";
matchedFlatFlammSeedScan::usage = "matchedFlatFlammSeedScan[n, seeds, k, radius, nbins] compares matched-flat and Flamm radial trends across seeds.";
matchedFlatFlammKScan::usage = "matchedFlatFlammKScan[n, seeds, ks, radius, nbins] repeats matched-flat/Flamm seed scans across k values.";
summarizeMatchedFlatKScan::usage = "summarizeMatchedFlatKScan[scan] summarizes matched-flat/Flamm correlations by k.";
runBHWithOriginalControl::usage = "runBHWithOriginalControl[label, f, targetK, seed, n, samplingMass, k, rg, nbins] runs a curved BH embedding with Flamm-style sampling plus the matched-flat control. Option \"UseBaseForCurved\" -> True forces the curved dataset to be the base Flamm graph; this is used for exact Schwarzschild-limit alignment.";
runRNOriginalControl::usage = "runRNOriginalControl[Q, seed, n, M, samplingMass, k, rg, nbins] runs one RN benchmark with the calibrated protocol.";
scanRNOriginalControl::usage = "scanRNOriginalControl[Qs, seeds, n, M, samplingMass, k, rg, nbins] runs RN across charges and seeds.";
scanBardeenOriginalControl::usage = "scanBardeenOriginalControl[gs, seeds, n, M, samplingMass, k, rg, nbins] runs Bardeen across parameters and seeds.";
scanHaywardOriginalControl::usage = "scanHaywardOriginalControl[ells, seeds, n, M, samplingMass, k, rg, nbins] runs Hayward across parameters and seeds.";
summarizeBHScan::usage = "summarizeBHScan[rows, parameterKey] summarizes BH-family seed scans.";

(* Export helpers *)
ensureDirectory::usage = "ensureDirectory[path] creates a directory if it does not exist.";
exportDatasetCSV::usage = "exportDatasetCSV[path, data] exports Dataset/list-of-associations data as CSV.";

Begin["`Private`"];

(* ================================================================ *)
(* Basic numerical helpers                                           *)
(* ================================================================ *)

ClearAll[validRealNumberQ, safeStd, safeCorrelation, rankData, spearmanCorr, cubicMeanDeviation];

validRealNumberQ[x_] :=
 Module[{y = Quiet[N[x]]},
  NumberQ[y] &&
   FreeQ[y, _Complex | Indeterminate | ComplexInfinity | DirectedInfinity[_]]
 ];

safeStd[x_List] :=
 Module[{vals = Select[N[x], validRealNumberQ]},
  If[Length[vals] < 2, 0, StandardDeviation[vals]]
 ];

safeCorrelation[x_List, y_List] :=
 Module[{pairs, xx, yy},
  pairs = Select[Transpose[{x, y}], validRealNumberQ[#[[1]]] && validRealNumberQ[#[[2]]] &];
  If[Length[pairs] < 3, Return[Missing["InsufficientData"]]];
  xx = N[pairs[[All, 1]]];
  yy = N[pairs[[All, 2]]];
  If[safeStd[xx] == 0 || safeStd[yy] == 0, Return[Missing["ZeroVariance"]]];
  Correlation[xx, yy]
 ];

rankData[x_List] :=
 Module[{pairs, groups, ranks, pos = 1, inds, r},
  pairs = SortBy[Transpose[{N[x], Range[Length[x]]}], First];
  groups = SplitBy[pairs, First];
  ranks = ConstantArray[0., Length[x]];
  Do[
   inds = group[[All, 2]];
   r = Mean[Range[pos, pos + Length[group] - 1]];
   ranks[[inds]] = r;
   pos += Length[group],
   {group, groups}
  ];
  ranks
 ];

spearmanCorr[x_List, y_List] :=
 Module[{pairs},
  pairs = Select[Transpose[{x, y}], validRealNumberQ[#[[1]]] && validRealNumberQ[#[[2]]] &];
  If[Length[pairs] < 3, Return[Missing["InsufficientData"]]];
  safeCorrelation[rankData[pairs[[All, 1]]], rankData[pairs[[All, 2]]]]
 ];

cubicMeanDeviation[x_List] :=
 Module[{vals = Select[N[x], validRealNumberQ], mu},
  If[Length[vals] < 2, Return[0.]];
  mu = Mean[vals];
  N[(Mean[Abs[vals - mu]^3])^(1/3)]
 ];

(* ================================================================ *)
(* Graph construction                                                *)
(* ================================================================ *)

ClearAll[distanceMatrixRadius, makeGeometricGraph, makeKNNGraph];

distanceMatrixRadius[pts_?MatrixQ, k_: 8, factor_: 1.15] :=
 Module[{dm, kdist},
  dm = DistanceMatrix[pts];
  kdist = Table[
    Sort[Delete[dm[[i]], i]][[Min[k, Length[pts] - 1]]],
    {i, Length[pts]}
   ];
  factor Median[kdist]
 ];

makeGeometricGraph[pts_?MatrixQ, eps_?NumericQ] :=
 Module[{n, dm, edges, coords, dim},
  n = Length[pts];
  dm = DistanceMatrix[pts];
  edges = Flatten[
    Table[
     If[i < j && dm[[i, j]] <= eps, UndirectedEdge[i, j], Nothing],
     {i, n}, {j, i + 1, n}
    ],
    1
   ];
  dim = Min[3, Length[pts[[1]]]];
  coords = Thread[Range[n] -> pts[[All, 1 ;; dim]]];
  Graph[
   Range[n],
   edges,
   VertexCoordinates -> coords,
   GraphLayout -> {"Dimension" -> dim},
   VertexSize -> Tiny,
   EdgeStyle -> Directive[Gray, Opacity[0.25]]
  ]
 ];

Options[makeKNNGraph] = {"VertexCoordinates" -> Automatic};

makeKNNGraph[points_?MatrixQ, k_Integer?Positive, OptionsPattern[]] :=
 Module[{n, nearest, edges, coords},
  n = Length[points];
  nearest = Nearest[points -> Range[n]];
  edges = DeleteDuplicates[
    Flatten[
     Table[(UndirectedEdge @@ Sort[{i, #}]) & /@ DeleteCases[nearest[points[[i]], k + 1], i], {i, n}],
     1
    ]
   ];
  coords = OptionValue["VertexCoordinates"];
  If[coords === Automatic,
   coords = If[Length[points[[1]]] >= 2, Thread[Range[n] -> points[[All, 1 ;; Min[3, Length[points[[1]]]]]]], Automatic]
  ];
  Graph[Range[n], edges, VertexCoordinates -> coords, VertexSize -> Tiny, EdgeStyle -> Directive[Gray, Opacity[0.35]]]
 ];

(* ================================================================ *)
(* Point clouds and dataset builders                                *)
(* ================================================================ *)

ClearAll[flatDiskPoints, hyperboloidPoints, flammPoints];

Options[flatDiskPoints] = {"RMin" -> 0.15, "RMax" -> 5.0};

flatDiskPoints[n_Integer?Positive, OptionsPattern[]] :=
 Module[{rmin, rmax, r, theta},
  rmin = OptionValue["RMin"];
  rmax = OptionValue["RMax"];
  r = Sqrt[RandomReal[{rmin^2, rmax^2}, n]];
  theta = RandomReal[{0, 2 Pi}, n];
  Transpose[{r Cos[theta], r Sin[theta]}]
 ];

Options[hyperboloidPoints] = {"RMin" -> 0.15, "RMax" -> 5.0, "HeightScale" -> 1.0};

hyperboloidPoints[n_Integer?Positive, OptionsPattern[]] :=
 Module[{rmin, rmax, scale, r, theta, z},
  rmin = OptionValue["RMin"];
  rmax = OptionValue["RMax"];
  scale = OptionValue["HeightScale"];
  r = Sqrt[RandomReal[{rmin^2, rmax^2}, n]];
  theta = RandomReal[{0, 2 Pi}, n];
  z = scale Sqrt[1 + r^2];
  Transpose[{r Cos[theta], r Sin[theta], z}]
 ];

Options[flammPoints] = {"RMin" -> Automatic, "RMax" -> 5.0};

flammPoints[n_Integer?Positive, M_: 1/2, OptionsPattern[]] :=
 Module[{rmin, rmax, r, theta, z},
  rmin = OptionValue["RMin"];
  rmax = OptionValue["RMax"];
  If[rmin === Automatic, rmin = 2 M + 0.2];
  r = Sqrt[RandomReal[{rmin^2, rmax^2}, n]];
  theta = RandomReal[{0, 2 Pi}, n];
  z = 2 Sqrt[2 M (r - 2 M)];
  Transpose[{r Cos[theta], r Sin[theta], z}]
 ];

ClearAll[schwarzschildK, buildFlatDataset, buildHyperboloidDataset, buildFlammDataset];

schwarzschildK[r_, M_: 1/2] := 48 M^2/r^6;

Options[buildFlatDataset] = {"RMin" -> 0.15, "RMax" -> 5.0};

buildFlatDataset[n_Integer?Positive, k_: 8, OptionsPattern[]] :=
 Module[{pts, r, eps, g},
  pts = flatDiskPoints[n, "RMin" -> OptionValue["RMin"], "RMax" -> OptionValue["RMax"]];
  r = Norm /@ pts;
  eps = distanceMatrixRadius[pts, k];
  g = makeGeometricGraph[pts, eps];
  <|
   "Type" -> "Flat",
   "Label" -> "Flat control",
   "N" -> n,
   "k" -> k,
   "Points" -> pts,
   "r" -> r,
   "Graph" -> g,
   "TargetK" -> ConstantArray[0., n],
   "Epsilon" -> eps
  |>
 ];

Options[buildHyperboloidDataset] = {"RMin" -> 0.15, "RMax" -> 5.0, "HeightScale" -> 1.0};

buildHyperboloidDataset[n_Integer?Positive, k_: 8, OptionsPattern[]] :=
 Module[{pts, r, eps, g},
  pts = hyperboloidPoints[
    n,
    "RMin" -> OptionValue["RMin"],
    "RMax" -> OptionValue["RMax"],
    "HeightScale" -> OptionValue["HeightScale"]
   ];
  r = Norm /@ pts[[All, 1 ;; 2]];
  eps = distanceMatrixRadius[pts, k];
  g = makeGeometricGraph[pts, eps];
  <|
   "Type" -> "Hyperboloid",
   "Label" -> "Hyperboloid-like benchmark",
   "N" -> n,
   "k" -> k,
   "Points" -> pts,
   "r" -> r,
   "Graph" -> g,
   "TargetK" -> ConstantArray[0., n],
   "Epsilon" -> eps
  |>
 ];

Options[buildFlammDataset] = {"RMin" -> Automatic, "RMax" -> 5.0};

buildFlammDataset[n_Integer?Positive, M_: 1/2, k_: 8, OptionsPattern[]] :=
 Module[{pts, r, targetK, eps, g},
  pts = flammPoints[n, M, "RMin" -> OptionValue["RMin"], "RMax" -> OptionValue["RMax"]];
  r = Norm /@ pts[[All, 1 ;; 2]];
  targetK = schwarzschildK[#, M] & /@ r;
  eps = distanceMatrixRadius[pts, k];
  g = makeGeometricGraph[pts, eps];
  <|
   "Type" -> "FlammSchwarzschild",
   "Label" -> "Flamm / Schwarzschild",
   "N" -> n,
   "k" -> k,
   "Mass" -> M,
   "M" -> M,
   "Points" -> pts,
   "r" -> r,
   "TargetK" -> targetK,
   "Graph" -> g,
   "Epsilon" -> eps
  |>
 ];

ClearAll[getKeySafe, getPointsSafe, getRadialCoordinates, getTargetKValues];

getKeySafe[data_, key_String, default_: Missing["KeyAbsent", key]] :=
 Module[{a, sym = ToExpression[key]},
  a = Which[
    AssociationQ[data], data,
    Head[data] === Dataset, Association[Normal[data]],
    True, Quiet[Check[Association[Normal[data]], <||>]]
   ];
  Which[
   AssociationQ[a] && KeyExistsQ[a, key], a[key],
   AssociationQ[a] && KeyExistsQ[a, sym], a[sym],
   True, default
  ]
 ];

getPointsSafe[data_] := getKeySafe[data, "Points", $Failed];
getRadialCoordinates[data_] := getKeySafe[data, "r", $Failed];
getTargetKValues[data_] := getKeySafe[data, "TargetK", $Failed];

buildMatchedFlatFromFlamm3[flammData_, k_Integer?Positive] :=
 Module[{pts3, xy, r, theta, pts, g},
  pts3 = getPointsSafe[flammData];
  If[pts3 === $Failed, Return[$Failed]];
  xy = pts3[[All, 1 ;; 2]];
  r = Norm /@ xy;
  theta = ArcTan @@@ xy;
  pts = Transpose[{r Cos[theta], r Sin[theta]}];
  g = makeKNNGraph[pts, k];
  <|
   "Type" -> "MatchedFlat",
   "Label" -> "Matched flat control (kNN)",
   "N" -> Length[pts],
   "k" -> k,
   "Points" -> pts,
   "r" -> r,
   "Graph" -> g,
   "TargetK" -> ConstantArray[0., Length[pts]]
  |>
 ];

buildMatchedFlatGeometricFromFlamm[flammData_, k_: 8] :=
 Module[{pts3, xy, r, theta, pts, eps, g},
  pts3 = getPointsSafe[flammData];
  If[pts3 === $Failed, Return[$Failed]];
  xy = pts3[[All, 1 ;; 2]];
  r = Norm /@ xy;
  theta = ArcTan @@@ xy;
  pts = Transpose[{r Cos[theta], r Sin[theta]}];
  eps = distanceMatrixRadius[pts, k];
  g = makeGeometricGraph[pts, eps];
  <|
   "Type" -> "MatchedFlatGeometric",
   "Label" -> "Matched flat control (epsilon graph)",
   "N" -> Length[pts],
   "k" -> k,
   "Points" -> pts,
   "r" -> r,
   "Graph" -> g,
   "TargetK" -> ConstantArray[0., Length[pts]],
   "Epsilon" -> eps
  |>
 ];

(* ================================================================ *)
(* Static spherical black-hole embeddings                            *)
(* ================================================================ *)

ClearAll[rnF, rnKretschmann, bardeenMass, bardeenF, haywardMass, haywardF];

rnF[M_, Q_][r_] := 1 - 2 M/r + Q^2/r^2;

rnKretschmann[M_, Q_][r_] := 48 M^2/r^6 - 96 M Q^2/r^7 + 56 Q^4/r^8;

bardeenMass[M_, g_][r_] := M r^3/(r^2 + g^2)^(3/2);
bardeenF[M_, g_][r_] := 1 - 2 bardeenMass[M, g][r]/r;

haywardMass[M_, ell_][r_] := M r^3/(r^3 + 2 M ell^2);
haywardF[M_, ell_][r_] := 1 - 2 haywardMass[M, ell][r]/r;

ClearAll[metricKretschmannFromF, metricKretschmannFunctionFromF];

metricKretschmannFromF[fexpr_, r_Symbol] :=
 Simplify[D[fexpr, {r, 2}]^2 + (2 D[fexpr, r]/r)^2 + 4 (1 - fexpr)^2/r^4];

metricKretschmannFunctionFromF[f_, hFactor_: 10^-4] :=
 Function[{r},
  Module[{h, fr, fp, fpp},
   h = hFactor Max[1, Abs[r]];
   fr = f[r];
   fp = (f[r + h] - f[r - h])/(2 h);
   fpp = (f[r + h] - 2 fr + f[r - h])/h^2;
   N[fpp^2 + (2 fp/r)^2 + 4 (1 - fr)^2/r^4]
  ]
 ];

ClearAll[phiFromCoordinates, embeddingZValuesFromF, buildStaticSphericalFromExistingSampling, buildStaticSphericalFromFlammSampling];

phiFromCoordinates[pts_?MatrixQ] := Map[ArcTan[#[[1]], #[[2]]] &, pts];

embeddingZValuesFromF[f_, rs_List] :=
 Module[{ord, sortedR, zprime, zp, dz, zsorted, zvals},
  ord = Ordering[rs];
  sortedR = N[rs[[ord]]];
  zprime[x_?NumericQ] := Sqrt[Max[0, 1/f[x] - 1]];
  zp = zprime /@ sortedR;
  dz = Differences[sortedR] ((Most[zp] + Rest[zp])/2);
  zsorted = Prepend[Accumulate[dz], 0.];
  zvals = ConstantArray[0., Length[rs]];
  zvals[[ord]] = zsorted;
  zvals
 ];

buildStaticSphericalFromExistingSampling[baseData_, f_, targetK_, k_: 8, label_: "StaticSpherical"] :=
 Module[{basePts, rs, phis, zs, pts, eps, g, kvals},
  basePts = getPointsSafe[baseData];
  rs = getRadialCoordinates[baseData];
  If[basePts === $Failed || rs === $Failed, Return[$Failed]];
  phis = phiFromCoordinates[basePts];
  zs = embeddingZValuesFromF[f, rs];
  pts = Transpose[{rs Cos[phis], rs Sin[phis], zs}];
  eps = distanceMatrixRadius[pts, k];
  g = makeGeometricGraph[pts, eps];
  kvals = Quiet[targetK /@ rs];
  <|
   "Type" -> label,
   "Label" -> label,
   "N" -> Length[pts],
   "k" -> k,
   "Points" -> pts,
   "r" -> rs,
   "theta" -> phis,
   "Graph" -> g,
   "TargetK" -> kvals,
   "Epsilon" -> eps
  |>
 ];

buildStaticSphericalFromFlammSampling[n_Integer?Positive, samplingMass_: 1/2, f_, targetK_, k_: 8, label_: "StaticSpherical"] :=
 Module[{base},
  base = buildFlammDataset[n, samplingMass, k];
  buildStaticSphericalFromExistingSampling[base, f, targetK, k, label]
 ];

buildRNDataset[n_Integer?Positive, M_: 1/2, Q_: 0.4, k_: 16] :=
 buildStaticSphericalFromFlammSampling[n, M, rnF[M, Q], rnKretschmann[M, Q], k, "RN_Q" <> ToString[Q]];

buildBardeenDataset[n_Integer?Positive, M_: 1/2, gpar_: 0.2, k_: 16] :=
 buildStaticSphericalFromFlammSampling[n, M, bardeenF[M, gpar], metricKretschmannFunctionFromF[bardeenF[M, gpar]], k, "Bardeen_g" <> ToString[gpar]];

buildHaywardDataset[n_Integer?Positive, M_: 1/2, ell_: 0.2, k_: 16] :=
 buildStaticSphericalFromFlammSampling[n, M, haywardF[M, ell], metricKretschmannFunctionFromF[haywardF[M, ell]], k, "Hayward_l" <> ToString[ell]];

buildBlackHoleFigureDatasets[n_Integer?Positive: 1000, M_: 1/2, k_: 16, Q_: 0.4, gpar_: 0.2, ell_: 0.2] :=
 Module[{base},
  base = buildFlammDataset[n, M, k];
  <|
   "Flamm" -> base,
   "RN" -> If[TrueQ[Chop[N[Q]] == 0],
     Join[base, <|"Type" -> "RN_Q0", "Label" -> "RN_Q0 exact Flamm limit"|>],
     buildStaticSphericalFromExistingSampling[base, rnF[M, Q], rnKretschmann[M, Q], k, "RN_Q" <> ToString[Q]]
    ],
   "Bardeen" -> If[TrueQ[Chop[N[gpar]] == 0],
     Join[base, <|"Type" -> "Bardeen_g0", "Label" -> "Bardeen_g0 exact Flamm limit"|>],
     buildStaticSphericalFromExistingSampling[base, bardeenF[M, gpar], metricKretschmannFunctionFromF[bardeenF[M, gpar]], k, "Bardeen_g" <> ToString[gpar]]
    ],
   "Hayward" -> If[TrueQ[Chop[N[ell]] == 0],
     Join[base, <|"Type" -> "Hayward_l0", "Label" -> "Hayward_l0 exact Flamm limit"|>],
     buildStaticSphericalFromExistingSampling[base, haywardF[M, ell], metricKretschmannFunctionFromF[haywardF[M, ell]], k, "Hayward_l" <> ToString[ell]]
    ]
  |>
 ];

(* ================================================================ *)
(* Basic plotting utilities                                          *)
(* ================================================================ *)

ClearAll[plotDataset2D, plotDataset3D, shellAndBallVertices, plotGeometryShell2D, plotGeometryShell3D, coloredGeometryPlot];

Options[plotDataset2D] = {"HighlightedVertices" -> {}, "PointSize" -> 0.010, "HighlightPointSize" -> 0.020};

plotDataset2D[data_, label_: Automatic, OptionsPattern[]] :=
 Module[{pts, hl, normal, title},
  pts = getPointsSafe[data];
  If[pts === $Failed, Return[$Failed]];
  pts = pts[[All, 1 ;; 2]];
  hl = OptionValue["HighlightedVertices"];
  normal = Complement[Range[Length[pts]], hl];
  title = If[label === Automatic, getKeySafe[data, "Label", "Dataset"], label];
  Show[
   ListPlot[pts[[normal]], AspectRatio -> 1, PlotStyle -> Directive[GrayLevel[0.55], PointSize[OptionValue["PointSize"]]], Frame -> True, Axes -> False, ImageSize -> Medium, PlotLabel -> title],
   If[Length[hl] > 0, ListPlot[pts[[hl]], AspectRatio -> 1, PlotStyle -> Directive[Red, PointSize[OptionValue["HighlightPointSize"]]]], Graphics[{}]]
  ]
 ];

Options[plotDataset3D] = {"HighlightedVertices" -> {}, "PointSize" -> 0.010, "HighlightPointSize" -> 0.025};

plotDataset3D[data_, label_: Automatic, OptionsPattern[]] :=
 Module[{pts, hl, normal, title},
  pts = getPointsSafe[data];
  If[pts === $Failed, Return[$Failed]];
  If[Length[pts[[1]]] < 3, Return[plotDataset2D[data, label]]];
  hl = OptionValue["HighlightedVertices"];
  normal = Complement[Range[Length[pts]], hl];
  title = If[label === Automatic, getKeySafe[data, "Label", "Dataset"], label];
  Show[
   ListPointPlot3D[pts[[normal]], PlotStyle -> Directive[GrayLevel[0.55], PointSize[OptionValue["PointSize"]]], BoxRatios -> Automatic, Axes -> True, ImageSize -> Medium, PlotLabel -> title],
   If[Length[hl] > 0, ListPointPlot3D[pts[[hl]], PlotStyle -> Directive[Red, PointSize[OptionValue["HighlightPointSize"]]]], Graphics3D[{}]]
  ]
 ];

shellAndBallVertices[data_, center_Integer, radius_Integer?Positive] :=
 Module[{g, verts, dist},
  g = getKeySafe[data, "Graph", $Failed];
  If[g === $Failed, Return[$Failed]];
  verts = VertexList[g];
  dist = AssociationThread[verts -> (GraphDistance[g, center, #] & /@ verts)];
  <|
   "Center" -> {center},
   "Ball" -> Select[verts, validRealNumberQ[dist[#]] && dist[#] <= radius &],
   "Shell" -> Select[verts, validRealNumberQ[dist[#]] && dist[#] == radius &]
  |>
 ];

plotGeometryShell2D[data_, center_Integer, radius_Integer?Positive, label_: Automatic] :=
 Module[{pts, sets, centerPts, shellPts, ballInteriorPts, otherPts, title},
  pts = getPointsSafe[data];
  If[pts === $Failed, Return[$Failed]];
  pts = pts[[All, 1 ;; 2]];
  sets = shellAndBallVertices[data, center, radius];
  If[sets === $Failed, Return[$Failed]];
  centerPts = sets["Center"];
  shellPts = Complement[sets["Shell"], centerPts];
  ballInteriorPts = Complement[sets["Ball"], shellPts, centerPts];
  otherPts = Complement[Range[Length[pts]], sets["Ball"]];
  title = If[label === Automatic, getKeySafe[data, "Label", "Dataset"], label];
  Show[
   ListPlot[pts[[otherPts]], PlotStyle -> Directive[GrayLevel[0.75], PointSize[0.008]]],
   ListPlot[pts[[ballInteriorPts]], PlotStyle -> Directive[LightBlue, PointSize[0.012]]],
   ListPlot[pts[[shellPts]], PlotStyle -> Directive[Orange, PointSize[0.016]]],
   ListPlot[pts[[centerPts]], PlotStyle -> Directive[Red, PointSize[0.025]]],
   Frame -> True, Axes -> False, AspectRatio -> 1, ImageSize -> Medium, PlotLabel -> title
  ]
 ];

plotGeometryShell3D[data_, center_Integer, radius_Integer?Positive, label_: Automatic] :=
 Module[{pts, sets, centerPts, shellPts, ballInteriorPts, otherPts, title},
  pts = getPointsSafe[data];
  If[pts === $Failed, Return[$Failed]];
  If[Length[pts[[1]]] < 3, Return[plotGeometryShell2D[data, center, radius, label]]];
  sets = shellAndBallVertices[data, center, radius];
  If[sets === $Failed, Return[$Failed]];
  centerPts = sets["Center"];
  shellPts = Complement[sets["Shell"], centerPts];
  ballInteriorPts = Complement[sets["Ball"], shellPts, centerPts];
  otherPts = Complement[Range[Length[pts]], sets["Ball"]];
  title = If[label === Automatic, getKeySafe[data, "Label", "Dataset"], label];
  Show[
   ListPointPlot3D[pts[[otherPts]], PlotStyle -> Directive[GrayLevel[0.75], PointSize[0.008]]],
   ListPointPlot3D[pts[[ballInteriorPts]], PlotStyle -> Directive[LightBlue, PointSize[0.012]]],
   ListPointPlot3D[pts[[shellPts]], PlotStyle -> Directive[Orange, PointSize[0.016]]],
   ListPointPlot3D[pts[[centerPts]], PlotStyle -> Directive[Red, PointSize[0.030]]],
   BoxRatios -> Automatic, ImageSize -> Medium, PlotLabel -> title
  ]
 ];

coloredGeometryPlot[data_, center_Integer, radius_Integer?Positive, label_: Automatic] :=
 Module[{pts = getPointsSafe[data]},
  If[pts === $Failed, Return[$Failed]];
  If[Length[pts[[1]]] >= 3, plotGeometryShell3D[data, center, radius, label], plotGeometryShell2D[data, center, radius, label]]
 ];

(* ================================================================ *)
(* Shortest-path anisotropy estimator                               *)
(* ================================================================ *)

ClearAll[precomputeAdjacency, shortestPathCountsFromAdj, pathAnisotropyEstimator, evaluateDataset];

precomputeAdjacency[g_Graph] := AssociationThread[VertexList[g] -> (AdjacencyList[g, #] & /@ VertexList[g])];

shortestPathCountsFromAdj[adj_Association, source_] :=
 Module[{verts, dist, counts, queue, v, nbrs},
  verts = Keys[adj];
  dist = AssociationThread[verts -> ConstantArray[Infinity, Length[verts]]];
  counts = AssociationThread[verts -> ConstantArray[0, Length[verts]]];
  dist[source] = 0;
  counts[source] = 1;
  queue = {source};
  While[Length[queue] > 0,
   v = First[queue];
   queue = Rest[queue];
   nbrs = adj[v];
   Do[
    Which[
     dist[u] === Infinity,
      dist[u] = dist[v] + 1;
      counts[u] = counts[v];
      queue = Append[queue, u],
     dist[u] == dist[v] + 1,
      counts[u] = counts[u] + counts[v]
    ],
    {u, nbrs}
   ];
  ];
  <|"Distance" -> dist, "Counts" -> counts|>
 ];

pathAnisotropyEstimator[g_Graph, p_Integer, radius_Integer?Positive, adj_: Automatic] :=
 Module[{adjacency, bfs, dist, counts, shell, countVals, logVals},
  adjacency = If[adj === Automatic, precomputeAdjacency[g], adj];
  bfs = shortestPathCountsFromAdj[adjacency, p];
  dist = bfs["Distance"];
  counts = bfs["Counts"];
  shell = Select[Keys[dist], dist[#] == radius &];
  If[Length[shell] < 2,
   Return[<|"Vertex" -> p, "Radius" -> radius, "ShellSize" -> Length[shell], "MeanCount" -> Missing["InsufficientShell"], "CMD" -> Missing["InsufficientShell"], "LogCMD" -> Missing["InsufficientShell"]|>]
  ];
  countVals = N[Lookup[counts, shell]];
  logVals = Log[countVals];
  <|
   "Vertex" -> p,
   "Radius" -> radius,
   "ShellSize" -> Length[shell],
   "MeanCount" -> Mean[countVals],
   "CMD" -> cubicMeanDeviation[countVals],
   "LogCMD" -> cubicMeanDeviation[logVals]
  |>
 ];

Options[evaluateDataset] = {"Vertices" -> Automatic};

evaluateDataset[data_Association, radius_Integer?Positive, OptionsPattern[]] :=
 Module[{g, verts, adj, rvals, kvals, rows},
  g = getKeySafe[data, "Graph", $Failed];
  If[g === $Failed, Return[$Failed]];
  verts = OptionValue["Vertices"];
  If[verts === Automatic, verts = VertexList[g]];
  adj = precomputeAdjacency[g];
  rvals = getRadialCoordinates[data];
  kvals = getTargetKValues[data];
  rows = Table[
    Module[{row = pathAnisotropyEstimator[g, v, radius, adj], extra = <||>},
     If[ListQ[rvals] && 1 <= v <= Length[rvals], extra = Join[extra, <|"r" -> rvals[[v]]|>]];
     If[ListQ[kvals] && 1 <= v <= Length[kvals], extra = Join[extra, <|"TargetK" -> kvals[[v]]|>]];
     Join[row, extra]
    ],
    {v, verts}
   ];
  Dataset[rows]
 ];

(* ================================================================ *)
(* Row cleaning/statistics                                           *)
(* ================================================================ *)

ClearAll[rowAssoc, getRowValue, cleanRows, estimatorStats];

rowAssoc[row_] :=
 Module[{n = Normal[row]},
  Which[
   AssociationQ[row], row,
   AssociationQ[n], n,
   ListQ[n] && AllTrue[n, MatchQ[#, _Rule] &], Association[n],
   True, <||>
  ]
 ];

getRowValue[row_, key_String] :=
 Module[{a = rowAssoc[row], sym = ToExpression[key]},
  Which[
   KeyExistsQ[a, key], a[key],
   KeyExistsQ[a, sym], a[sym],
   True, Missing["KeyAbsent", key]
  ]
 ];

cleanRows[res_, estimator_: "LogCMD"] :=
 Module[{rows = Normal[res]},
  Reap[
    Do[
     Module[{v, rad, e, k, r, shell},
      v = getRowValue[row, "Vertex"];
      rad = getRowValue[row, "Radius"];
      e = Quiet[N[getRowValue[row, estimator]]];
      k = Quiet[N[getRowValue[row, "TargetK"]]];
      r = Quiet[N[getRowValue[row, "r"]]];
      shell = getRowValue[row, "ShellSize"];
      If[validRealNumberQ[e] && validRealNumberQ[k],
       Sow[<|"Vertex" -> v, "Radius" -> rad, "ShellSize" -> shell, "r" -> r, "Estimator" -> N[e], "TargetK" -> N[k]|>]
      ]
     ],
     {row, rows}
    ]
   ][[2]] /. {} -> {{}} // First
 ];

estimatorStats[res_, estimator_: "LogCMD"] :=
 Module[{rows, e, k},
  rows = cleanRows[res, estimator];
  If[Length[rows] < 3, Return[<|"Rows" -> Length[rows], "Pearson" -> Missing["InsufficientData"], "Spearman" -> Missing["InsufficientData"]|>]];
  e = Lookup[rows, "Estimator"];
  k = Lookup[rows, "TargetK"];
  <|"Rows" -> Length[rows], "Pearson" -> safeCorrelation[e, k], "Spearman" -> spearmanCorr[e, k]|>
 ];

(* ================================================================ *)
(* Estimator-colored 3D plots                                        *)
(* ================================================================ *)

ClearAll[estimatorValuesByVertex, allGoodEstimatorValues, redBlueCF, blueRedCF, estimatorColoredGraph3DSharedScale, makeBlackHoleFamilyPanel];

estimatorValuesByVertex[data_, res_, estimator_: "LogCMD"] :=
 Module[{vals, rows, v, e},
  vals = ConstantArray[Missing["NoValue"], Length[getPointsSafe[data]]];
  rows = Normal[res];
  Do[
   v = getRowValue[row, "Vertex"];
   e = Quiet[N[getRowValue[row, estimator]]];
   If[IntegerQ[v] && 1 <= v <= Length[vals] && NumericQ[e], vals[[v]] = e],
   {row, rows}
  ];
  vals
 ];

allGoodEstimatorValues[pairs_] := Flatten[Table[Cases[estimatorValuesByVertex[pair[[1]], pair[[2]], "LogCMD"], _?NumericQ], {pair, pairs}]];

redBlueCF = Function[t, Blend[{RGBColor[0.75, 0.05, 0.05], RGBColor[0.95, 0.95, 0.98], RGBColor[0.08, 0.22, 0.92]}, t]];
blueRedCF = Function[t, Blend[{RGBColor[0.08, 0.22, 0.92], RGBColor[0.95, 0.95, 0.98], RGBColor[0.75, 0.05, 0.05]}, t]];

Options[estimatorColoredGraph3DSharedScale] = {"ColorFunction" -> redBlueCF, "MaxEdges" -> All};

estimatorColoredGraph3DSharedScale[data_, res_, range_, title_: "", pointSize_: 0.012, edgeOpacity_: 0.10, OptionsPattern[]] :=
 Module[{pts, g, vals, cf, color, edges, usedEdges, edgePrimitives, pointPrimitives},
  pts = getPointsSafe[data];
  g = getKeySafe[data, "Graph", $Failed];
  vals = estimatorValuesByVertex[data, res, "LogCMD"];
  cf = OptionValue["ColorFunction"];
  color[v_] := If[NumericQ[v], cf[Clip[Rescale[v, range], {0, 1}]], GrayLevel[0.75]];
  edges = EdgeList[g];
  usedEdges = If[OptionValue["MaxEdges"] === All || Length[edges] <= OptionValue["MaxEdges"], edges, RandomSample[edges, OptionValue["MaxEdges"]]];
  edgePrimitives = {Directive[GrayLevel[0.25], Opacity[edgeOpacity], Thin], Map[Function[e, With[{i = First[List @@ e], j = Last[List @@ e]}, Line[{pts[[i]], pts[[j]]}]]], usedEdges]};
  pointPrimitives = MapThread[{color[#2], Point[#1]} &, {pts, vals}];
  Graphics3D[{edgePrimitives, PointSize[pointSize], pointPrimitives}, Boxed -> False, Axes -> False, PlotRange -> All, ImageSize -> 420, PlotLabel -> Style[title, 15], ViewPoint -> {2.2, -2.4, 1.3}, Background -> White]
 ];

makeBlackHoleFamilyPanel[pairs_, labels_, range_: Automatic, colorFunction_: redBlueCF] :=
 Module[{vals, sharedRange, graphics},
  vals = allGoodEstimatorValues[pairs];
  sharedRange = If[range === Automatic, MinMax[vals], range];
  graphics = MapThread[estimatorColoredGraph3DSharedScale[#1[[1]], #1[[2]], sharedRange, #2, 0.012, 0.10, "ColorFunction" -> colorFunction] &, {pairs, labels}];
  Legended[GraphicsGrid[Partition[graphics, 2], Spacings -> {0.2, 0.4}, ImageSize -> Large], BarLegend[{colorFunction, sharedRange}, LegendLabel -> "C_log"]]
 ];

(* ================================================================ *)
(* Radial binning and scans                                          *)
(* ================================================================ *)

ClearAll[assignRadialBins, radialBinComparisonNoClean, radialBinComparisonFast2, binnedFlammComparison, corrFromBinned, curvatureScoreBins, canonicalCurvatureBinnedCorrelations];

assignRadialBins[rows_List, nbins_Integer?Positive] :=
 Module[{rmin, rmax},
  If[Length[rows] == 0, Return[{}]];
  rmin = Min[Lookup[rows, "r"]];
  rmax = Max[Lookup[rows, "r"]];
  If[rmax == rmin, Return[Map[Join[#, <|"Bin" -> 1|>] &, rows]]];
  Map[Function[row, Join[row, <|"Bin" -> Min[nbins, Max[1, 1 + Floor[nbins (row["r"] - rmin)/(rmax - rmin)]]]|>]], rows]
 ];

radialBinComparisonNoClean[data_, res_, estimator_: "LogCMD", nbins_Integer?Positive] :=
 Module[{rawRows, rvals, rows, withBins, grouped, bins},
  rawRows = Normal[res];
  rvals = getRadialCoordinates[data];
  If[rvals === $Failed, Return[$Failed]];
  rows = Reap[
      Do[
       Module[{v, e},
        v = getRowValue[row, "Vertex"];
        e = Quiet[N[getRowValue[row, estimator]]];
        If[IntegerQ[v] && 1 <= v <= Length[rvals] && validRealNumberQ[e], Sow[<|"Vertex" -> v, "r" -> N[rvals[[v]]], "Estimator" -> N[e]|>]]
       ],
       {row, rawRows}
      ]
     ][[2]] /. {} -> {{}} // First;
  If[Length[rows] == 0, Return[{}]];
  withBins = assignRadialBins[rows, nbins];
  grouped = GroupBy[withBins, #"Bin" &];
  bins = Sort[Keys[grouped]];
  Table[<|"Bin" -> b, "rMean" -> Mean[Lookup[grouped[b], "r"]], "EstimatorMean" -> Mean[Lookup[grouped[b], "Estimator"]], "EstimatorStd" -> safeStd[Lookup[grouped[b], "Estimator"]], "Count" -> Length[grouped[b]]|>, {b, bins}]
 ];

radialBinComparisonFast2[data_, res_, estimator_: "LogCMD", nbins_Integer?Positive] :=
 Module[{rows, withBins, grouped, bins},
  rows = cleanRows[res, estimator];
  rows = Select[rows, validRealNumberQ[#"r"] && validRealNumberQ[#"Estimator"] &];
  If[Length[rows] == 0, Return[{}]];
  withBins = assignRadialBins[rows, nbins];
  grouped = GroupBy[withBins, #"Bin" &];
  bins = Sort[Keys[grouped]];
  Table[<|"Bin" -> b, "rMean" -> Mean[Lookup[grouped[b], "r"]], "EstimatorMean" -> Mean[Lookup[grouped[b], "Estimator"]], "EstimatorStd" -> safeStd[Lookup[grouped[b], "Estimator"]], "Count" -> Length[grouped[b]]|>, {b, bins}]
 ];

binnedFlammComparison[data_, res_, estimator_: "LogCMD", nbins_Integer?Positive] :=
 Module[{rows, withBins, grouped, bins},
  rows = cleanRows[res, estimator];
  rows = Select[rows, validRealNumberQ[#"r"] && validRealNumberQ[#"Estimator"] && validRealNumberQ[#"TargetK"] &];
  If[Length[rows] == 0, Return[{}]];
  withBins = assignRadialBins[rows, nbins];
  grouped = GroupBy[withBins, #"Bin" &];
  bins = Sort[Keys[grouped]];
  Table[<|"Bin" -> b, "rMean" -> Mean[Lookup[grouped[b], "r"]], "KMean" -> Mean[Lookup[grouped[b], "TargetK"]], "EstimatorMean" -> Mean[Lookup[grouped[b], "Estimator"]], "EstimatorStd" -> safeStd[Lookup[grouped[b], "Estimator"]], "Count" -> Length[grouped[b]]|>, {b, bins}]
 ];

corrFromBinned[binned_List] := safeCorrelation[Lookup[binned, "rMean"], Lookup[binned, "EstimatorMean"]];

(*
  Reporting-level curvature-oriented convention.

  LogCMD remains the raw estimator C_log.  For the black-hole-family
  summary table, it is clearer to report

      C_curv = - C_log

  so that larger C_curv corresponds to a stronger inward curvature signal
  under the fixed graph-construction protocol.  This does not modify the
  estimator; it only adds a signed reporting column to binned summaries.
*)

curvatureScoreBins[binned_List] :=
 Map[
  Function[row,
   If[
    AssociationQ[row] && KeyExistsQ[row, "EstimatorMean"] && validRealNumberQ[row["EstimatorMean"]],
    Join[row, <|"CurvatureScoreMean" -> -N[row["EstimatorMean"]]|>],
    row
    ]
   ],
  binned
  ];

canonicalCurvatureBinnedCorrelations[binned_List] :=
 Module[{scoredBins, goodBins, rvals, cvals, kvals, corrR, corrK},

  scoredBins = curvatureScoreBins[binned];

  goodBins =
   SortBy[
    Select[
     scoredBins,
     AssociationQ[#] &&
       KeyExistsQ[#, "rMean"] &&
       KeyExistsQ[#, "CurvatureScoreMean"] &&
       validRealNumberQ[#["rMean"]] &&
       validRealNumberQ[#["CurvatureScoreMean"]] &
     ],
    #["rMean"] &
    ];

  If[Length[goodBins] < 3,
   Return[
    <|
     "Convention" -> "Ccurv = -C_log; rMean increases outward",
     "BinsUsed" -> Length[goodBins],
     "CorrR_Ccurv" -> Missing["InsufficientData"],
     "CorrLogK_Ccurv" -> Missing["InsufficientData"]
     |>
    ]
   ];

  rvals = N @ Lookup[goodBins, "rMean"];
  cvals = N @ Lookup[goodBins, "CurvatureScoreMean"];
  corrR = N @ safeCorrelation[rvals, cvals];

  corrK =
   If[
    AllTrue[
     goodBins,
     KeyExistsQ[#, "KMean"] &&
       validRealNumberQ[#["KMean"]] &&
       #["KMean"] > 0 &
     ],
    kvals = N @ Lookup[goodBins, "KMean"];
    N @ safeCorrelation[Log[kvals], cvals],
    Missing["NoPositiveKMean"]
    ];

  <|
   "Convention" -> "Ccurv = -C_log; rMean increases outward",
   "BinsUsed" -> Length[goodBins],
   "CorrR_Ccurv" -> corrR,
   "CorrLogK_Ccurv" -> corrK
   |>
  ];

ClearAll[rowsWithRAndEstimator, commonRadialBinnedComparison, plotCommonRadialProfiles];

rowsWithRAndEstimator[data_, res_, estimator_: "LogCMD"] :=
 Module[{rawRows, rvals, rows},
  rawRows = Normal[res];
  rvals = getRadialCoordinates[data];
  If[rvals === $Failed, Return[{}]];
  rows = Reap[
      Do[
       Module[{v = getRowValue[row, "Vertex"], e = Quiet[N[getRowValue[row, estimator]]]},
        If[IntegerQ[v] && 1 <= v <= Length[rvals] && validRealNumberQ[e], Sow[<|"Vertex" -> v, "r" -> N[rvals[[v]]], "Estimator" -> N[e]|>]]
       ],
       {row, rawRows}
      ]
     ][[2]] /. {} -> {{}} // First;
  rows
 ];

commonRadialBinnedComparison[data1_, res1_, data2_, res2_, estimator_: "LogCMD", nbins_Integer?Positive: 12] :=
 Module[{rows1, rows2, rmin, rmax, assignBin, addBins, grouped1, grouped2, commonBins},
  rows1 = rowsWithRAndEstimator[data1, res1, estimator];
  rows2 = rowsWithRAndEstimator[data2, res2, estimator];
  If[Length[rows1] == 0 || Length[rows2] == 0, Return[{}]];
  rmin = Max[Min[Lookup[rows1, "r"]], Min[Lookup[rows2, "r"]]];
  rmax = Min[Max[Lookup[rows1, "r"]], Max[Lookup[rows2, "r"]]];
  assignBin[r_] := Min[nbins, Max[1, 1 + Floor[nbins (r - rmin)/(rmax - rmin)]]];
  addBins[rows_] := Select[Map[Function[row, Join[row, <|"Bin" -> assignBin[row["r"]]|>]], rows], rmin <= #"r" <= rmax &];
  grouped1 = GroupBy[addBins[rows1], #"Bin" &];
  grouped2 = GroupBy[addBins[rows2], #"Bin" &];
  commonBins = Intersection[Keys[grouped1], Keys[grouped2]];
  N[Table[<|"Bin" -> b, "rMean" -> Mean[Join[Lookup[grouped1[b], "r"], Lookup[grouped2[b], "r"]]], "EstimatorMean1" -> Mean[Lookup[grouped1[b], "Estimator"]], "EstimatorMean2" -> Mean[Lookup[grouped2[b], "Estimator"]], "Count1" -> Length[grouped1[b]], "Count2" -> Length[grouped2[b]]|>, {b, Sort[commonBins]}]]
 ];

plotCommonRadialProfiles[commonBinned_, label1_: "Matched flat", label2_: "Flamm", title_: "Common-bin radial path anisotropy", file_: None] :=
 Module[{plot},
  plot = ListLinePlot[
    {Transpose[{Lookup[commonBinned, "rMean"], Lookup[commonBinned, "EstimatorMean1"]}], Transpose[{Lookup[commonBinned, "rMean"], Lookup[commonBinned, "EstimatorMean2"]}]},
    Frame -> True,
    FrameLabel -> {"Common radial bin coordinate r", "Binned mean LogCMD"},
    PlotLegends -> {label1, label2},
    PlotMarkers -> Automatic,
    GridLines -> Automatic,
    ImageSize -> Large,
    PlotLabel -> title
   ];
  If[file =!= None, Export[file, plot]];
  plot
 ];

(* ================================================================ *)
(* Scans                                                            *)
(* ================================================================ *)

ClearAll[radiusScan, fastRadiusScan, multiSeedRadiusScan, summarizeSeedScan];

radiusScan[data_Association, radii_List, estimator_: "LogCMD"] :=
 Dataset[Table[Module[{res, stats}, res = evaluateDataset[data, radius]; stats = estimatorStats[res, estimator]; Join[<|"Type" -> getKeySafe[data, "Type", Missing["Type"]], "N" -> getKeySafe[data, "N", Missing["N"]], "k" -> getKeySafe[data, "k", Missing["k"]], "Radius" -> radius|>, stats]], {radius, radii}]];
fastRadiusScan[data_Association, radii_List, estimator_: "LogCMD"] := radiusScan[data, radii, estimator];

multiSeedRadiusScan[n_Integer?Positive, seeds_List, radii_List, k_Integer?Positive, M_: 1/2, estimator_: "LogCMD"] :=
 Dataset[Flatten[Table[Module[{data, scan}, SeedRandom[seed]; data = buildFlammDataset[n, M, k]; scan = Normal[radiusScan[data, radii, estimator]]; Map[Join[<|"Seed" -> seed|>, #] &, scan]], {seed, seeds}], 1]];

summarizeSeedScan[scan_] :=
 Module[{rows, groups},
  rows = Normal[scan];
  groups = GroupBy[rows, {#"N" &, #"k" &, #"Radius" &}];
  Dataset[KeyValueMap[Function[{key, vals}, <|"N" -> key[[1]], "k" -> key[[2]], "Radius" -> key[[3]], "MeanPearson" -> Mean[Select[Lookup[vals, "Pearson"], validRealNumberQ]], "MeanSpearman" -> Mean[Select[Lookup[vals, "Spearman"], validRealNumberQ]], "StdSpearman" -> safeStd[Lookup[vals, "Spearman"]], "MeanRows" -> Mean[Lookup[vals, "Rows"]], "Seeds" -> Length[vals]|>], groups]]
 ];

ClearAll[matchedFlatFlammSeedScan, matchedFlatFlammKScan, summarizeMatchedFlatKScan];

matchedFlatFlammSeedScan[n_Integer?Positive, seeds_List, k_Integer?Positive, radius_Integer?Positive, nbins_Integer?Positive, M_: 1/2] :=
 Table[
  Module[{flamm, matchedFlat, resFlamm, resMatched, binnedFlamm, binnedMatched, corrFlamm, corrMatched},
   SeedRandom[seed];
   flamm = buildFlammDataset[n, M, k];
   matchedFlat = buildMatchedFlatFromFlamm3[flamm, k];
   resFlamm = evaluateDataset[flamm, radius];
   resMatched = evaluateDataset[matchedFlat, radius];
   binnedFlamm = radialBinComparisonNoClean[flamm, resFlamm, "LogCMD", nbins];
   binnedMatched = radialBinComparisonNoClean[matchedFlat, resMatched, "LogCMD", nbins];
   corrFlamm = corrFromBinned[binnedFlamm];
   corrMatched = corrFromBinned[binnedMatched];
   <|"Seed" -> seed, "N" -> n, "k" -> k, "Radius" -> radius, "Bins" -> nbins, "MatchedFlatCorrRLogCMD" -> N[corrMatched], "FlammCorrRLogCMD" -> N[corrFlamm], "Difference" -> N[corrFlamm - corrMatched]|>
  ],
  {seed, seeds}
 ];

matchedFlatFlammKScan[n_Integer?Positive, seeds_List, ks_List, radius_Integer?Positive, nbins_Integer?Positive, M_: 1/2] := Flatten[Table[matchedFlatFlammSeedScan[n, seeds, k, radius, nbins, M], {k, ks}], 1];

summarizeMatchedFlatKScan[scan_List] :=
 Module[{ks = Sort[DeleteDuplicates[Lookup[scan, "k"]]]},
  Table[Module[{rows, mf, fl}, rows = Select[scan, #"k" == k &]; mf = Lookup[rows, "MatchedFlatCorrRLogCMD"]; fl = Lookup[rows, "FlammCorrRLogCMD"]; <|"k" -> k, "MatchedFlatMean" -> Mean[mf], "MatchedFlatStd" -> safeStd[mf], "FlammMean" -> Mean[fl], "FlammStd" -> safeStd[fl]|>], {k, ks}]
 ];

ClearAll[runBHWithOriginalControl, runRNOriginalControl, scanRNOriginalControl, scanBardeenOriginalControl, scanHaywardOriginalControl, summarizeBHScan];

Options[runBHWithOriginalControl] = {"UseBaseForCurved" -> False};

runBHWithOriginalControl[label_, f_, targetK_, seed_: 1234, n_: 1000, samplingMass_: 1/2, k_: 16, rg_: 3, nbins_: 12, OptionsPattern[]] :=
 Module[{base, curved, flat, resCurved, resFlat, bCurved, bFlat, bCurvedCcurv, bFlatCcurv, statsCurved, statsFlat, corrCurved, corrFlat, corrLogK},
  SeedRandom[seed];
  base = buildFlammDataset[n, samplingMass, k];

  (*
    Critical alignment convention:
    The black-hole-family scans are calibrated against the original
    Schwarzschild/Flamm construction. For zero-deformation limits
    (RN Q=0, Bardeen g=0, Hayward ell=0), the curved dataset is set
    equal to the base Flamm graph itself. This prevents accidental
    mismatch from numerical re-embedding or graph-construction changes.
  *)
  curved =
   If[TrueQ[OptionValue["UseBaseForCurved"]],
    Join[base, <|"Type" -> label, "Label" -> label|>],
    buildStaticSphericalFromExistingSampling[base, f, targetK, k, label]
   ];

  (* The matched-flat control is always built from the same base Flamm sampling. *)
  flat = buildMatchedFlatFromFlamm3[base, k];
  resCurved = evaluateDataset[curved, rg];
  resFlat = evaluateDataset[flat, rg];
  bCurved = binnedFlammComparison[curved, resCurved, "LogCMD", nbins];
  bFlat = radialBinComparisonNoClean[flat, resFlat, "LogCMD", nbins];

  (*
    Journal/reporting convention:
      C_curv = - C_log.
    The raw LogCMD estimator is not changed.  We only use the signed
    CurvatureScoreMean column for the BH-family summary correlations.
  *)
  bCurvedCcurv = curvatureScoreBins[bCurved];
  bFlatCcurv = curvatureScoreBins[bFlat];
  statsCurved = canonicalCurvatureBinnedCorrelations[bCurved];
  statsFlat = canonicalCurvatureBinnedCorrelations[bFlat];
  corrCurved = statsCurved["CorrR_Ccurv"];
  corrFlat = statsFlat["CorrR_Ccurv"];
  corrLogK = statsCurved["CorrLogK_Ccurv"];

  <|
   "Label" -> label,
   "Seed" -> seed,
   "N" -> n,
   "k" -> k,
   "Radius" -> rg,
   "Bins" -> nbins,
   "CcurvDefinition" -> "Ccurv = -C_log",
   "Curved" -> curved,
   "Flat" -> flat,
   "CurvedResult" -> resCurved,
   "FlatResult" -> resFlat,
   "CurvedBinned" -> bCurved,
   "FlatBinned" -> bFlat,
   "CurvedBinnedCcurv" -> bCurvedCcurv,
   "FlatBinnedCcurv" -> bFlatCcurv,
   "CorrR_Ccurv" -> N[corrCurved],
   "CorrLogK_Ccurv" -> N[corrLogK],
   "MatchedFlatCorrR_Ccurv" -> N[corrFlat],
   "Difference_Ccurv" -> N[corrCurved - corrFlat]
   |>
 ];

runRNOriginalControl[Q_, seed_: 1234, n_: 1000, M_: 1/2, samplingMass_: 1/2, k_: 16, rg_: 3, nbins_: 12] :=
 Join[
  <|"Q" -> Q|>,
  runBHWithOriginalControl[
   "RN_Q" <> ToString[Q],
   rnF[M, Q],
   rnKretschmann[M, Q],
   seed, n, samplingMass, k, rg, nbins,
   "UseBaseForCurved" -> TrueQ[Chop[N[Q]] == 0]
  ]
 ];

scanRNOriginalControl[Qs_List, seeds_List, n_: 1000, M_: 1/2, samplingMass_: 1/2, k_: 16, rg_: 3, nbins_: 12] :=
 Flatten[Table[KeyDrop[runRNOriginalControl[Q, seed, n, M, samplingMass, k, rg, nbins], {"Curved", "Flat", "CurvedResult", "FlatResult", "CurvedBinned", "FlatBinned", "CurvedBinnedCcurv", "FlatBinnedCcurv"}], {Q, Qs}, {seed, seeds}], 1];

scanBardeenOriginalControl[gs_List, seeds_List, n_: 1000, M_: 1/2, samplingMass_: 1/2, k_: 16, rg_: 3, nbins_: 12] :=
 Flatten[
  Table[
   Module[
    {out = runBHWithOriginalControl[
       "Bardeen_g" <> ToString[g],
       bardeenF[M, g],
       metricKretschmannFunctionFromF[bardeenF[M, g]],
       seed, n, samplingMass, k, rg, nbins,
       "UseBaseForCurved" -> TrueQ[Chop[N[g]] == 0]
     ]},
    Join[<|"g" -> g|>, KeyDrop[out, {"Curved", "Flat", "CurvedResult", "FlatResult", "CurvedBinned", "FlatBinned", "CurvedBinnedCcurv", "FlatBinnedCcurv"}]]
   ],
   {g, gs}, {seed, seeds}
  ],
  1
 ];

scanHaywardOriginalControl[ells_List, seeds_List, n_: 1000, M_: 1/2, samplingMass_: 1/2, k_: 16, rg_: 3, nbins_: 12] :=
 Flatten[
  Table[
   Module[
    {out = runBHWithOriginalControl[
       "Hayward_l" <> ToString[ell],
       haywardF[M, ell],
       metricKretschmannFunctionFromF[haywardF[M, ell]],
       seed, n, samplingMass, k, rg, nbins,
       "UseBaseForCurved" -> TrueQ[Chop[N[ell]] == 0]
     ]},
    Join[<|"ell" -> ell|>, KeyDrop[out, {"Curved", "Flat", "CurvedResult", "FlatResult", "CurvedBinned", "FlatBinned", "CurvedBinnedCcurv", "FlatBinnedCcurv"}]]
   ],
   {ell, ells}, {seed, seeds}
  ],
  1
 ];

summarizeBHScan[rows_List, parameterKey_String] :=
 Module[{groups},
  groups = GroupBy[rows, #[parameterKey] &];
  Dataset[
   KeyValueMap[
    Function[{par, vals},
     <|
      parameterKey -> par,
      "CcurvDefinition" -> "Ccurv = -C_log",
      "CorrR_CcurvMean" -> Mean[Select[Lookup[vals, "CorrR_Ccurv"], validRealNumberQ]],
      "CorrR_CcurvStd" -> safeStd[Lookup[vals, "CorrR_Ccurv"]],
      "CorrLogK_CcurvMean" -> Mean[Select[Lookup[vals, "CorrLogK_Ccurv"], validRealNumberQ]],
      "CorrLogK_CcurvStd" -> safeStd[Lookup[vals, "CorrLogK_Ccurv"]],
      "MatchedFlatCorrR_CcurvMean" -> Mean[Select[Lookup[vals, "MatchedFlatCorrR_Ccurv"], validRealNumberQ]],
      "MatchedFlatCorrR_CcurvStd" -> safeStd[Lookup[vals, "MatchedFlatCorrR_Ccurv"]],
      "Difference_CcurvMean" -> Mean[Select[Lookup[vals, "Difference_Ccurv"], validRealNumberQ]],
      "Difference_CcurvStd" -> safeStd[Lookup[vals, "Difference_Ccurv"]],
      "Seeds" -> Length[vals]
      |>
     ],
    groups
    ]
   ]
 ];

(* ================================================================ *)
(* Export helpers                                                    *)
(* ================================================================ *)

ClearAll[ensureDirectory, exportDatasetCSV];

ensureDirectory[path_String] := If[! DirectoryQ[path], CreateDirectory[path, CreateIntermediateDirectories -> True]];

exportDatasetCSV[path_String, data_] :=
 Module[{dir = DirectoryName[path], rows = Normal[data]},
  If[StringQ[dir] && dir =!= "", ensureDirectory[dir]];
  Export[path, rows]
 ];

End[];
EndPackage[];
