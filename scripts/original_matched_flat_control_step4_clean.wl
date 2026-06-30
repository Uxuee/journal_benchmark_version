(* ::Package:: *)
(* Original matched-flat control for black-hole graph diagnostics. *)

(* ::Section:: *)
(* 0. Clean settings *)

ClearAll["Global`*"];

(* ----------------------------- *)
(* Main run settings             *)
(* ----------------------------- *)

M = 1.0;

(* Current-scale version used in the new Step 4 work. *)
rMin = 2.2;
rMax = 12.0;

(* Old-paper scale, if you want to reproduce the old Flamm-only experiment exactly. *)
(* M = 1/2; rMin = 1.05; rMax = 5.0; *)

nSamples = 1000;
randomSeed = 1234;

kNN = 16;
epsilonFactor = 1.15;
fixedGraphShellRadius = 3;
nRadialBins = 12;

(* Use All for the real run. Use 5 or 10 only for a smoke test. *)
maxRefsPerBin = All;
(* maxRefsPerBin = 5; *)

qRN = 0.5;
gBardeen = 0.5;
lHayward = 0.5;

runTag = "originalMatchedFlatControl_randomN" <> ToString[nSamples] <>
  "_k" <> ToString[kNN] <>
  "_epsFactor" <> StringReplace[ToString[epsilonFactor], "." -> "p"] <>
  "_rg" <> ToString[fixedGraphShellRadius];

baseDir = If[$Notebooks, NotebookDirectory[], Directory[]];

outputDir = FileNameJoin[{baseDir, "step4_outputs_" <> runTag}];
checkpointDir = FileNameJoin[{baseDir, "step4_checkpoints_" <> runTag}];

If[! DirectoryQ[outputDir], CreateDirectory[outputDir]];
If[! DirectoryQ[checkpointDir], CreateDirectory[checkpointDir]];

runTag

(* ::Section:: *)
(* 1. Metric functions *)

ClearAll[
  fSchwarzschild,
  fRN,
  mBardeen,
  fBardeen,
  mHayward,
  fHayward
];

fSchwarzschild[r_?NumericQ] := 1 - 2 M/r;

fRN[r_?NumericQ] := 1 - 2 M/r + qRN^2/r^2;

mBardeen[r_?NumericQ] := M r^3/(r^2 + gBardeen^2)^(3/2);
fBardeen[r_?NumericQ] := 1 - 2 mBardeen[r]/r;

mHayward[r_?NumericQ] := M r^3/(r^3 + 2 M lHayward^2);
fHayward[r_?NumericQ] := 1 - 2 mHayward[r]/r;

(* ::Section:: *)
(* 2. Embedding and sampling utilities *)

ClearAll[
  safeZPrime,
  makeZInterpolation,
  makeRandomRPhiList,
  makeGridRPhiList
];

safeZPrime[f_][r_?NumericQ] := Module[{fr},
  fr = f[r];
  If[NumericQ[fr] && fr > 0,
    Sqrt[Max[0, 1/fr - 1]],
    0
  ]
];

makeZInterpolation[f_, rmin_, rmax_, nGrid_: 3000] := Module[
  {rs, vals, increments, zs},
  rs = N @ Subdivide[rmin, rmax, nGrid];
  vals = safeZPrime[f] /@ rs;
  increments = Differences[rs] * MovingAverage[vals, 2];
  zs = Prepend[Accumulate[increments], 0];
  Interpolation[Transpose[{rs, zs}], InterpolationOrder -> 1]
];

makeRandomRPhiList[n_Integer, rmin_, rmax_, seed_Integer: 1234] := Module[
  {rs, phis},
  SeedRandom[seed];
  rs = RandomReal[{rmin, rmax}, n];
  phis = RandomReal[{0, 2 Pi}, n];
  Transpose[{rs, phis}]
];

makeGridRPhiList[nR_Integer, nPhi_Integer, rmin_, rmax_] := Module[
  {rs, phis},
  rs = N @ Subdivide[rmin, rmax, nR - 1];
  phis = N @ Most[Subdivide[0, 2 Pi, nPhi]];
  Flatten[Table[{r, phi}, {r, rs}, {phi, phis}], 1]
];

(* ::Section:: *)
(* 3. Original graph-construction protocol *)

ClearAll[
  distanceMatrixRadius,
  makeGeometricGraph,
  makeKNNGraph,
  makeSurfaceDataFromRPhi,
  makeFlatMatchedDataOriginal
];

(* Original radius rule: epsilon = factor * Median[k-th neighbor distance]. *)
distanceMatrixRadius[pts_List, k_Integer: 8, factor_: 1.15] := Module[
  {dm, sortedRows, kthDistances},
  dm = DistanceMatrix[pts];
  sortedRows = Sort /@ dm;
  (* First distance is zero, the point itself; k+1 gives the k-th other neighbor. *)
  kthDistances = sortedRows[[All, Min[k + 1, Length[pts]]]];
  N[factor * Median[kthDistances]]
];

makeGeometricGraph[pts_List, eps_?NumericQ] := Module[
  {n, dm, edges},
  n = Length[pts];
  dm = DistanceMatrix[pts];
  edges = Flatten[
    Table[
      If[i < j && dm[[i, j]] <= eps,
        UndirectedEdge[i, j],
        Nothing
      ],
      {i, n}, {j, n}
    ],
    1
  ];
  SimpleGraph[Graph[Range[n], edges]]
];

makeKNNGraph[pts_List, k_Integer] := Module[
  {n, nf, pairs, sortedPairs, edges},
  n = Length[pts];
  nf = Nearest[pts -> Range[n]];
  pairs = Flatten[
    Table[
      Thread[i -> DeleteCases[nf[pts[[i]], k + 1], i]],
      {i, n}
    ],
    1
  ];
  sortedPairs = Sort /@ (List @@@ pairs);
  edges = DeleteDuplicates[UndirectedEdge @@@ sortedPairs];
  SimpleGraph[Graph[Range[n], edges]]
];

makeSurfaceDataFromRPhi[name_String, f_, rPhiList_List] := Module[
  {zfun, pts, vertices, radiusAssociation, phiAssociation,
   coordinateAssociation, eps, graph},

  zfun = makeZInterpolation[f, rMin, rMax];

  pts = ({#[[1]] Cos[#[[2]]], #[[1]] Sin[#[[2]]], zfun[#[[1]]]} &) /@ rPhiList;

  vertices = Range[Length[pts]];

  radiusAssociation = AssociationThread[vertices -> N[rPhiList[[All, 1]]]];
  phiAssociation = AssociationThread[vertices -> N[rPhiList[[All, 2]]]];
  coordinateAssociation = AssociationThread[vertices -> pts];

  eps = distanceMatrixRadius[pts, kNN, epsilonFactor];
  graph = makeGeometricGraph[pts, eps];

  <|
    "Name" -> name,
    "Graph" -> graph,
    "Points" -> pts,
    "RadiusAssociation" -> radiusAssociation,
    "PhiAssociation" -> phiAssociation,
    "CoordinateAssociation" -> coordinateAssociation,
    "GraphConstruction" -> "radiusThresholdOriginal",
    "GraphK" -> kNN,
    "GraphEpsilon" -> eps,
    "EpsilonFactor" -> epsilonFactor,
    "Sampling" -> "sharedRandomPolar"
  |>
];

(* Original matched-flat convention: project same sample to z=0 and use direct kNN. *)
makeFlatMatchedDataOriginal[curvedData_Association] := Module[
  {name, vertices, radiusAssociation, phiAssociation,
   radii, phis, pts, coordinateAssociation, graph},

  name = curvedData["Name"] <> "_matchedFlatOriginal";

  radiusAssociation = curvedData["RadiusAssociation"];
  phiAssociation = curvedData["PhiAssociation"];
  vertices = Keys[radiusAssociation];

  radii = Lookup[radiusAssociation, vertices];
  phis = Lookup[phiAssociation, vertices];

  pts = MapThread[
    {#1 Cos[#2], #1 Sin[#2], 0.0} &,
    {radii, phis}
  ];

  coordinateAssociation = AssociationThread[vertices -> pts];
  graph = makeKNNGraph[pts, kNN];

  <|
    "Name" -> name,
    "Graph" -> graph,
    "Points" -> pts,
    "RadiusAssociation" -> radiusAssociation,
    "PhiAssociation" -> phiAssociation,
    "CoordinateAssociation" -> coordinateAssociation,
    "GraphConstruction" -> "matchedFlatOriginalKNN",
    "GraphK" -> kNN,
    "GraphEpsilon" -> Missing["NotUsed"],
    "EpsilonFactor" -> Missing["NotUsed"],
    "Sampling" -> "sameRphiAsCurved"
  |>
];

(* ::Section:: *)
(* 4. Build curved graphs and original matched-flat controls *)

baseRPhiList = makeRandomRPhiList[nSamples, rMin, rMax, randomSeed];

curvedGraphData = {
  makeSurfaceDataFromRPhi["SchwarzschildFlamm", fSchwarzschild, baseRPhiList],
  makeSurfaceDataFromRPhi["ReissnerNordstrom", fRN, baseRPhiList],
  makeSurfaceDataFromRPhi["Bardeen", fBardeen, baseRPhiList],
  makeSurfaceDataFromRPhi["Hayward", fHayward, baseRPhiList]
};

flatGraphData = makeFlatMatchedDataOriginal /@ curvedGraphData;

allGraphData = Join[curvedGraphData, flatGraphData];

graphSummary = <|
    "Name" -> #["Name"],
    "Construction" -> #["GraphConstruction"],
    "k" -> #["GraphK"],
    "Epsilon" -> #["GraphEpsilon"],
    "EpsilonFactor" -> #["EpsilonFactor"],
    "Vertices" -> VertexCount[#["Graph"]],
    "Edges" -> EdgeCount[#["Graph"]],
    "ConnectedQ" -> ConnectedGraphQ[#["Graph"]]
  |> & /@ allGraphData;

Export[FileNameJoin[{outputDir, "graph_summary_original_control.csv"}], graphSummary];
Export[FileNameJoin[{outputDir, "graph_summary_original_control.mx"}], graphSummary];

Dataset[graphSummary]

(* ::Section:: *)
(* 5. Radial bins and reference vertices *)

ClearAll[
  makeRadialBinAssociation,
  selectReferenceVertices
];

makeRadialBinAssociation[radiusAssociation_Association, nBins_Integer] := Module[
  {vertices, radii, rlo, rhi},
  vertices = Keys[radiusAssociation];
  radii = Values[radiusAssociation];
  rlo = Min[radii];
  rhi = Max[radii];
  AssociationThread[
    vertices ->
      (
        Clip[
          1 + Floor[nBins (Lookup[radiusAssociation, #] - rlo)/(rhi - rlo)],
          {1, nBins}
        ] & /@ vertices
      )
  ]
];

selectReferenceVertices[radiusAssociation_Association, nBins_Integer, maxPerBin_: All] := Module[
  {binAssociation, vertices, grouped},
  vertices = Keys[radiusAssociation];
  binAssociation = makeRadialBinAssociation[radiusAssociation, nBins];
  grouped = GroupBy[vertices, binAssociation[#] &];
  If[maxPerBin === All,
    Flatten[Values[grouped]],
    SeedRandom[1234];
    Flatten[RandomSample[#, Min[Length[#], maxPerBin]] & /@ Values[grouped]]
  ]
];

(* ::Section:: *)
(* 6. Shortest-path multiplicity diagnostics *)

ClearAll[
  localClusteringCoefficient,
  shortestPathCountsFromSource,
  sourceDiagnostics,
  diagnosticsForGraph
];

localClusteringCoefficient[g_Graph, v_] := Module[
  {neighbors, d, neighborPairs, existingEdges},
  neighbors = AdjacencyList[g, v];
  d = Length[neighbors];
  If[d < 2, Return[0.]];
  neighborPairs = Subsets[neighbors, {2}];
  existingEdges = Count[neighborPairs, pair_ /; EdgeQ[g, UndirectedEdge @@ pair]];
  N[2 existingEdges/(d (d - 1))]
];

shortestPathCountsFromSource[g_Graph, source_, adjacencyAssociation_Association] := Module[
  {vertices, distances, reachable, ordered, sigma, predecessors},

  vertices = VertexList[g];

  distances = AssociationThread[
    vertices,
    GraphDistance[g, source, #] & /@ vertices
  ];

  reachable = Select[
    vertices,
    NumericQ[distances[#]] && distances[#] < Infinity &
  ];

  ordered = SortBy[reachable, distances[#] &];

  sigma = AssociationThread[vertices, ConstantArray[0, Length[vertices]]];
  sigma[source] = 1;

  Do[
    If[v =!= source,
      predecessors = Select[
        Lookup[adjacencyAssociation, v, {}],
        KeyExistsQ[distances, #] &&
          NumericQ[distances[#]] &&
          distances[#] == distances[v] - 1 &
      ];
      sigma[v] = Total[Lookup[sigma, predecessors, 0]];
    ],
    {v, ordered}
  ];

  sigma
];

sourceDiagnostics[
   g_Graph,
   source_,
   adjacencyAssociation_Association,
   rg_Integer : fixedGraphShellRadius
   ] := Module[
  {vertices, distances, reachable, shellVertices, shellSize,
   sigma, logSigma, cLog, meanLogMultiplicity, stdLogMultiplicity,
   degree, clustering},

  vertices = DeleteCases[VertexList[g], source];

  distances = AssociationThread[
    vertices,
    GraphDistance[g, source, #] & /@ vertices
  ];

  reachable = Select[
    vertices,
    NumericQ[distances[#]] && distances[#] < Infinity &
  ];

  shellVertices = Select[reachable, distances[#] == rg &];
  shellSize = Length[shellVertices];

  sigma = shortestPathCountsFromSource[g, source, adjacencyAssociation];

  logSigma = If[shellSize > 0,
    Log[N[Lookup[sigma, shellVertices]]],
    {}
  ];

  cLog = If[Length[logSigma] > 1,
    Mean[Abs[logSigma - Mean[logSigma]]^3]^(1/3),
    0.
  ];

  meanLogMultiplicity = If[Length[logSigma] > 0,
    Mean[logSigma],
    Missing["EmptyShell"]
  ];

  stdLogMultiplicity = If[Length[logSigma] > 1,
    StandardDeviation[logSigma],
    0.
  ];

  degree = VertexDegree[g, source];
  clustering = localClusteringCoefficient[g, source];

  <|
   "GraphShellRadius" -> N[rg],
   "CLog" -> N[cLog],
   "MeanShellSize" -> N[shellSize],
   "StdShellSize" -> 0.,
   "MeanLogShellSize" -> N[Log[1 + shellSize]],
   "StdLogShellSize" -> 0.,
   "MeanLogMultiplicityByShell" -> N[meanLogMultiplicity],
   "StdLogMultiplicityByShell" -> N[stdLogMultiplicity],
   "Degree" -> N[degree],
   "ClusteringCoefficient" -> N[clustering]
   |>
  ];

diagnosticsForGraph[data_Association, nBins_Integer, maxPerBin_: All] := Module[
  {g, name, radiusAssociation, binAssociation, refs,
   adjacencyAssociation, rows},

  name = data["Name"];
  g = data["Graph"];
  radiusAssociation = data["RadiusAssociation"];

  binAssociation = makeRadialBinAssociation[radiusAssociation, nBins];
  refs = selectReferenceVertices[radiusAssociation, nBins, maxPerBin];

  adjacencyAssociation = AssociationThread[
    VertexList[g],
    AdjacencyList[g, #] & /@ VertexList[g]
  ];

  rows = Table[
    Join[
      <|
        "GraphName" -> name,
        "GraphConstruction" -> data["GraphConstruction"],
        "SourceVertex" -> v,
        "Radius" -> N[radiusAssociation[v]],
        "RadialBin" -> N[binAssociation[v]]
      |>,
      sourceDiagnostics[g, v, adjacencyAssociation, fixedGraphShellRadius]
    ],
    {v, refs}
  ];

  rows
];

(* ::Section:: *)
(* 7. Aggregation helpers *)

ClearAll[
  numericValues,
  meanOrMissing,
  sdOrZero,
  aggregateRowsByGraphAndBin,
  numericizeRows
];

numericValues[x_List] := Select[N[x], NumericQ];

meanOrMissing[x_List] := Module[{vals},
  vals = numericValues[x];
  If[Length[vals] > 0, N[Mean[vals]], Missing["NoNumericData"]]
];

sdOrZero[x_List] := Module[{vals},
  vals = numericValues[x];
  If[Length[vals] > 1, N[StandardDeviation[vals]], 0.]
];

aggregateRowsByGraphAndBin[rows_List] := Module[
  {groups, diagnosticKeys},

  diagnosticKeys = {
    "CLog",
    "MeanShellSize",
    "StdShellSize",
    "MeanLogShellSize",
    "StdLogShellSize",
    "MeanLogMultiplicityByShell",
    "StdLogMultiplicityByShell",
    "Degree",
    "ClusteringCoefficient"
  };

  groups = Values @ GroupBy[rows, {#["GraphName"], #["RadialBin"]} &] /. Key[x_] :> x;

  Table[
    Association[
      Flatten[
        {
          "GraphName" -> First[group]["GraphName"],
          "GraphConstruction" -> First[group]["GraphConstruction"],
          "RadialBin" -> First[group]["RadialBin"],
          "RadiusMean" -> N[Mean[Lookup[group, "Radius"]]],
          "NReferenceVertices" -> Length[group],

          Table[
            {
              key <> "_Mean" -> meanOrMissing[Lookup[group, key]],
              key <> "_SD" -> sdOrZero[Lookup[group, key]]
            },
            {key, diagnosticKeys}
          ]
        }
      ]
    ],
    {group, groups}
  ]
];

numericizeRows[rows_List] := (Association @ KeyValueMap[#1 -> Quiet@Check[N[#2], #2] &, #]) & /@ rows;

(* ::Section:: *)
(* 8. Checkpointed execution *)

ClearAll[runOneGraphCheckpoint];

runOneGraphCheckpoint[data_Association, nBins_Integer, maxPerBin_] := Module[
  {name, refTag, sourceFile, binnedFile, time, rows, binnedRows},

  name = data["Name"];
  refTag = If[maxPerBin === All, "refsAll", "refsPerBin" <> ToString[maxPerBin]];

  sourceFile = FileNameJoin[{checkpointDir, name <> "_" <> refTag <> "_sourceRows.mx"}];
  binnedFile = FileNameJoin[{checkpointDir, name <> "_" <> refTag <> "_binnedRows.mx"}];

  If[FileExistsQ[sourceFile] && FileExistsQ[binnedFile],
    Print["Already exists, loading: ", name];
    Return[
      <|
        "Name" -> name,
        "SourceRows" -> Import[sourceFile],
        "BinnedRows" -> Import[binnedFile]
      |>
    ];
  ];

  Print["Running diagnostics for: ", name];

  {time, rows} = AbsoluteTiming[
    CheckAbort[
      diagnosticsForGraph[data, nBins, maxPerBin],
      $Aborted
    ]
  ];

  If[rows === $Aborted || rows === $Failed,
    Print["FAILED or aborted: ", name];
    Return[$Failed];
  ];

  rows = numericizeRows[rows];

  Print["Finished ", name, " in ", NumberForm[time/60, {5, 2}], " minutes."];

  Export[sourceFile, rows];

  binnedRows = numericizeRows[aggregateRowsByGraphAndBin[rows]];

  Export[binnedFile, binnedRows];

  <|
    "Name" -> name,
    "SourceRows" -> rows,
    "BinnedRows" -> binnedRows
  |>
];

(* ::Section:: *)
(* 9. Smoke test on Schwarzschild/Flamm *)

smokeMaxRefsPerBin = 5;

smokeFlamm = runOneGraphCheckpoint[
  SelectFirst[allGraphData, #["Name"] == "SchwarzschildFlamm" &],
  nRadialBins,
  smokeMaxRefsPerBin
];

Dataset[smokeFlamm["BinnedRows"]]

(* ::Section:: *)
(* 10. Full run *)

(* Run this after the smoke test looks reasonable. With maxRefsPerBin = All it may take time. *)

resultsByGraph = Table[
  runOneGraphCheckpoint[data, nRadialBins, maxRefsPerBin],
  {data, allGraphData}
];

refTag = If[maxRefsPerBin === All, "refsAll", "refsPerBin" <> ToString[maxRefsPerBin]];
sourceFiles = FileNames["*_" <> refTag <> "_sourceRows.mx", checkpointDir];
binnedFiles = FileNames["*_" <> refTag <> "_binnedRows.mx", checkpointDir];

allSourceRows = Join @@ (Import /@ sourceFiles);
allBinnedRows = Join @@ (Import /@ binnedFiles);

allSourceRows = numericizeRows[allSourceRows];
allBinnedRows = numericizeRows[allBinnedRows];

Export[FileNameJoin[{outputDir, "source_vertex_diagnostics_original_control.mx"}], allSourceRows];
Export[FileNameJoin[{outputDir, "radial_bin_diagnostics_original_control.mx"}], allBinnedRows];
Export[FileNameJoin[{outputDir, "radial_bin_diagnostics_original_control.csv"}], allBinnedRows];

Dataset[allBinnedRows]

(* ::Section:: *)
(* 11. Profiles, plots, and Schwarzschild/control gap tables *)

ClearAll[
  diagnosticsToCompare,
  getBinnedProfile,
  plotBHvsFlat,
  correlationSafe,
  corrWithRadius,
  baselineGapTable,
  meanAbsoluteProfileGap,
  standardizedProfileGap,
  profileShapeCorrelation
];

diagnosticsToCompare = {
  "CLog",
  "MeanLogShellSize",
  "MeanLogMultiplicityByShell",
  "Degree",
  "ClusteringCoefficient"
};

getBinnedProfile[graphName_String, diagnosticKey_String] := Module[
  {rows, key},
  key = diagnosticKey <> "_Mean";
  rows = Select[allBinnedRows, #["GraphName"] == graphName &];
  SortBy[
    ({#["RadiusMean"], #[key]} &) /@ rows,
    First
  ]
];

plotBHvsFlat[bhName_String, flatName_String, diagnosticKey_String] := Module[
  {bh, flat, plt, outFile},
  bh = getBinnedProfile[bhName, diagnosticKey];
  flat = getBinnedProfile[flatName, diagnosticKey];
  plt = ListLinePlot[
    {bh, flat},
    Frame -> True,
    FrameLabel -> {"r", diagnosticKey},
    PlotLegends -> {bhName, flatName},
    PlotMarkers -> Automatic,
    ImageSize -> Large,
    PlotLabel -> diagnosticKey <> ": black-hole graph vs original matched-flat control"
  ];
  outFile = FileNameJoin[{outputDir, "profile_" <> diagnosticKey <> "_" <> bhName <> "_vs_originalFlat.pdf"}];
  Export[outFile, plt];
  plt
];

correlationSafe[x_List, y_List] := Module[{pairs},
  pairs = Select[Transpose[{x, y}], VectorQ[#, NumericQ] &];
  If[Length[pairs] >= 3,
    N[Correlation[pairs[[All, 1]], pairs[[All, 2]]]],
    Missing["TooFewPoints"]
  ]
];

corrWithRadius[graphName_String, diagnosticKey_String] := Module[{profile},
  profile = getBinnedProfile[graphName, diagnosticKey];
  correlationSafe[profile[[All, 1]], profile[[All, 2]]]
];

baselineGapTable[bhName_String, flatName_String, diagnosticKeys_List] := Table[
  With[
    {
      bhCorr = corrWithRadius[bhName, key],
      flatCorr = corrWithRadius[flatName, key]
    },
    <|
      "Geometry" -> bhName,
      "Control" -> flatName,
      "Diagnostic" -> key,
      "BHCorrR" -> bhCorr,
      "FlatCorrR" -> flatCorr,
      "Gap" -> N[bhCorr - flatCorr],
      "AbsGap" -> N[Abs[bhCorr - flatCorr]]
    |>
  ],
  {key, diagnosticKeys}
];

meanAbsoluteProfileGap[bhName_String, flatName_String, diagnosticKey_String] := Module[
  {bh, flat, commonR, bhAssoc, flatAssoc, diffs},
  bh = getBinnedProfile[bhName, diagnosticKey];
  flat = getBinnedProfile[flatName, diagnosticKey];
  bhAssoc = Association[Rule @@@ bh];
  flatAssoc = Association[Rule @@@ flat];
  commonR = Intersection[Keys[bhAssoc], Keys[flatAssoc]];
  diffs = Abs[N[Lookup[bhAssoc, commonR]] - N[Lookup[flatAssoc, commonR]]];
  If[Length[diffs] > 0, N[Mean[diffs]], Missing["NoCommonBins"]]
];

standardizedProfileGap[bhName_String, flatName_String, diagnosticKey_String] := Module[
  {bh, flat, commonR, bhAssoc, flatAssoc, bhVals, flatVals, pooledSD},
  bh = getBinnedProfile[bhName, diagnosticKey];
  flat = getBinnedProfile[flatName, diagnosticKey];
  bhAssoc = Association[Rule @@@ bh];
  flatAssoc = Association[Rule @@@ flat];
  commonR = Intersection[Keys[bhAssoc], Keys[flatAssoc]];
  bhVals = N[Lookup[bhAssoc, commonR]];
  flatVals = N[Lookup[flatAssoc, commonR]];
  pooledSD = StandardDeviation[Join[bhVals, flatVals]];
  If[pooledSD > 0, N[Mean[Abs[bhVals - flatVals]]/pooledSD], Missing["ZeroPooledSD"]]
];

profileShapeCorrelation[bhName_String, flatName_String, diagnosticKey_String] := Module[
  {bh, flat, commonR, bhAssoc, flatAssoc, bhVals, flatVals},
  bh = getBinnedProfile[bhName, diagnosticKey];
  flat = getBinnedProfile[flatName, diagnosticKey];
  bhAssoc = Association[Rule @@@ bh];
  flatAssoc = Association[Rule @@@ flat];
  commonR = Intersection[Keys[bhAssoc], Keys[flatAssoc]];
  bhVals = N[Lookup[bhAssoc, commonR]];
  flatVals = N[Lookup[flatAssoc, commonR]];
  If[Length[bhVals] > 2, N[Correlation[bhVals, flatVals]], Missing["TooFewPoints"]]
];

flammFlatName = "SchwarzschildFlamm_matchedFlatOriginal";

Do[
  plotBHvsFlat["SchwarzschildFlamm", flammFlatName, key],
  {key, diagnosticsToCompare}
];

flammGapTable = baselineGapTable[
  "SchwarzschildFlamm",
  flammFlatName,
  diagnosticsToCompare
];

flammAmplitudeGapTable = Table[
  <|
    "Diagnostic" -> key,
    "MeanAbsProfileGap" -> meanAbsoluteProfileGap["SchwarzschildFlamm", flammFlatName, key],
    "StandardizedProfileGap" -> standardizedProfileGap["SchwarzschildFlamm", flammFlatName, key],
    "BHFlatShapeCorrelation" -> profileShapeCorrelation["SchwarzschildFlamm", flammFlatName, key]
  |>,
  {key, diagnosticsToCompare}
];

Export[FileNameJoin[{outputDir, "flamm_correlation_gap_table_original_control.csv"}], flammGapTable];
Export[FileNameJoin[{outputDir, "flamm_amplitude_shape_gap_table_original_control.csv"}], flammAmplitudeGapTable];

Dataset[flammGapTable]

(* ::Section:: *)
(* 12. Amplitude/shape gap table *)

Dataset[flammAmplitudeGapTable]

(* ::Section:: *)
(* 13. Optional: run all geometry-vs-original-flat gap tables *)

allGeometryGapTables = Flatten[
  Table[
    baselineGapTable[
      geom,
      geom <> "_matchedFlatOriginal",
      diagnosticsToCompare
    ],
    {geom, {"SchwarzschildFlamm", "ReissnerNordstrom", "Bardeen", "Hayward"}}
  ],
  1
];

allGeometryAmplitudeGapTables = Flatten[
  Table[
    Table[
      <|
        "Geometry" -> geom,
        "Control" -> geom <> "_matchedFlatOriginal",
        "Diagnostic" -> key,
        "MeanAbsProfileGap" -> meanAbsoluteProfileGap[geom, geom <> "_matchedFlatOriginal", key],
        "StandardizedProfileGap" -> standardizedProfileGap[geom, geom <> "_matchedFlatOriginal", key],
        "BHFlatShapeCorrelation" -> profileShapeCorrelation[geom, geom <> "_matchedFlatOriginal", key]
      |>,
      {key, diagnosticsToCompare}
    ],
    {geom, {"SchwarzschildFlamm", "ReissnerNordstrom", "Bardeen", "Hayward"}}
  ],
  1
];

Export[FileNameJoin[{outputDir, "all_geometry_correlation_gap_table_original_control.csv"}], allGeometryGapTables];
Export[FileNameJoin[{outputDir, "all_geometry_amplitude_shape_gap_table_original_control.csv"}], allGeometryAmplitudeGapTables];

Dataset[allGeometryGapTables]

(* ::Section:: *)
(* 14. Final sensitivity checks for manuscript *)

(*
  This section contains the final reproducibility checks used in the journal
  manuscript:
    - graph summary table,
    - graph-shell-radius sweep for CLog,
    - original matched-flat KNN vs strict radius-threshold matched-flat control.

  Use All for the final run. For a quick smoke test, set rgSweepMaxRefsPerBin
  or sensitivityMaxRefsPerBin to a small integer such as 20.
*)

extraDir = FileNameJoin[{outputDir, "extra_sensitivity_checks"}];
If[! DirectoryQ[extraDir], CreateDirectory[extraDir]];

ClearAll[graphSummaryExtended];

graphSummaryExtended[data_Association] := Module[
  {g = data["Graph"], deg},
  deg = VertexDegree[g];
  <|
    "Name" -> data["Name"],
    "Construction" -> data["GraphConstruction"],
    "k" -> Lookup[data, "GraphK", Missing["NotAvailable"]],
    "Epsilon" -> Lookup[data, "GraphEpsilon", Missing["NotUsed"]],
    "EpsilonFactor" -> Lookup[data, "EpsilonFactor", Missing["NotUsed"]],
    "Vertices" -> VertexCount[g],
    "Edges" -> EdgeCount[g],
    "ConnectedQ" -> ConnectedGraphQ[g],
    "MeanDegree" -> N@Mean[deg],
    "MinDegree" -> Min[deg],
    "MaxDegree" -> Max[deg],
    "EdgeDensity" -> N@GraphDensity[g]
  |>
];

graphSummaryTableExtended = graphSummaryExtended /@ allGraphData;
Export[
  FileNameJoin[{extraDir, "graph_summary_extended_original_control.csv"}],
  graphSummaryTableExtended
];

ClearAll[
  diagnosticsForGraphAtRG,
  getProfileFromBinnedRows,
  correlationSafeLocal,
  corrWithRadiusFromBinnedRows,
  gapTableFromBinnedRows
];

diagnosticsForGraphAtRG[
   data_Association,
   nBins_Integer : nRadialBins,
   maxPerBin_ : All,
   rg_Integer : fixedGraphShellRadius
] := Module[
  {g, name, radiusAssociation, binAssociation, refs, adjacencyAssociation, rows},

  name = data["Name"];
  g = data["Graph"];
  radiusAssociation = data["RadiusAssociation"];
  binAssociation = makeRadialBinAssociation[radiusAssociation, nBins];
  refs = selectReferenceVertices[radiusAssociation, nBins, maxPerBin];

  adjacencyAssociation = AssociationThread[
    VertexList[g],
    AdjacencyList[g, #] & /@ VertexList[g]
  ];

  rows = Table[
    Join[
      <|
        "GraphName" -> name,
        "GraphConstruction" -> data["GraphConstruction"],
        "SourceVertex" -> v,
        "Radius" -> N[radiusAssociation[v]],
        "RadialBin" -> N[binAssociation[v]],
        "GraphShellRadiusRequested" -> rg
      |>,
      sourceDiagnostics[g, v, adjacencyAssociation, rg]
    ],
    {v, refs}
  ];

  rows
];

getProfileFromBinnedRows[rows_, graphName_String, diagnostic_String] := Module[
  {sub, key},
  key = diagnostic <> "_Mean";
  sub = Select[
    rows,
    #["GraphName"] == graphName &&
      NumericQ[#["RadiusMean"]] &&
      NumericQ[Lookup[#, key, Missing[]]] &
  ];
  SortBy[({#["RadiusMean"], Lookup[#, key]} & /@ sub), First]
];

correlationSafeLocal[pairs_] := Module[{clean},
  clean = Select[pairs, VectorQ[#, NumericQ] &];
  If[Length[clean] < 3,
    Missing["TooFewPoints"],
    N@Correlation[clean[[All, 1]], clean[[All, 2]]]
  ]
];

corrWithRadiusFromBinnedRows[rows_, graphName_String, diagnostic_String] :=
  correlationSafeLocal[getProfileFromBinnedRows[rows, graphName, diagnostic]];

gapTableFromBinnedRows[
   rows_,
   bhName_String,
   flatName_String,
   diagnostics_List,
   rg_Integer : fixedGraphShellRadius
] := Table[
  With[
    {
      cBH = corrWithRadiusFromBinnedRows[rows, bhName, diag],
      cF = corrWithRadiusFromBinnedRows[rows, flatName, diag]
    },
    <|
      "Geometry" -> bhName,
      "Control" -> flatName,
      "Diagnostic" -> diag,
      "BHCorrR" -> cBH,
      "FlatCorrR" -> cF,
      "Gap" -> If[NumericQ[cBH] && NumericQ[cF], N[cBH - cF], Missing["NotAvailable"]],
      "AbsGap" -> If[NumericQ[cBH] && NumericQ[cF], N[Abs[cBH - cF]], Missing["NotAvailable"]],
      "GraphShellRadius" -> rg
    |>
  ],
  {diag, diagnostics}
];

(* Shell-radius sweep for Schwarzschild/Flamm under the original matched-flat control. *)

rgSweepValues = {2, 3, 4, 5};
rgSweepMaxRefsPerBin = All;
diagnosticsForRGSweep = diagnosticsToCompare;

bhSchwarzschild = First@Select[allGraphData, #["Name"] == "SchwarzschildFlamm" &];
flatSchwarzschildOriginal =
  First@Select[allGraphData, #["Name"] == "SchwarzschildFlamm_matchedFlatOriginal" &];

ClearAll[runRGSweepOne];
runRGSweepOne[rg_Integer] := Module[
  {rowsBH, rowsFlat, rowsAll, binned, gaps},
  Print["Running rg = ", rg];

  rowsBH = diagnosticsForGraphAtRG[bhSchwarzschild, nRadialBins, rgSweepMaxRefsPerBin, rg];
  rowsFlat = diagnosticsForGraphAtRG[flatSchwarzschildOriginal, nRadialBins, rgSweepMaxRefsPerBin, rg];

  rowsAll = Join[rowsBH, rowsFlat];
  binned = numericizeRows[aggregateRowsByGraphAndBin[rowsAll]];

  Export[
    FileNameJoin[{extraDir, "rg_sweep_binned_rows_rg" <> ToString[rg] <> ".csv"}],
    binned
  ];

  gaps = gapTableFromBinnedRows[
    binned,
    "SchwarzschildFlamm",
    "SchwarzschildFlamm_matchedFlatOriginal",
    diagnosticsForRGSweep,
    rg
  ];

  Export[
    FileNameJoin[{extraDir, "rg_sweep_gap_rows_rg" <> ToString[rg] <> ".csv"}],
    gaps
  ];

  gaps
];

rgSweepGapTable = Flatten[runRGSweepOne /@ rgSweepValues, 1];
Export[
  FileNameJoin[{extraDir, "rg_sweep_gap_table_schwarzschild_original_control.csv"}],
  rgSweepGapTable
];

rgSweepCLogOnly = Select[rgSweepGapTable, #["Diagnostic"] == "CLog" &];
Export[
  FileNameJoin[{extraDir, "rg_sweep_clog_only_schwarzschild_original_control.csv"}],
  rgSweepCLogOnly
];

(* Strict radius-threshold matched-flat control. *)

ClearAll[makeFlatMatchedDataRadiusThreshold];

makeFlatMatchedDataRadiusThreshold[curvedData_Association] := Module[
  {name, vertices, radiusAssociation, phiAssociation, radii, phis, pts,
   coordinateAssociation, eps, graph},

  name = curvedData["Name"] <> "_matchedFlatStrictRadius";
  radiusAssociation = curvedData["RadiusAssociation"];
  phiAssociation = curvedData["PhiAssociation"];
  vertices = Keys[radiusAssociation];

  radii = Lookup[radiusAssociation, vertices];
  phis = Lookup[phiAssociation, vertices];

  pts = MapThread[{#1 Cos[#2], #1 Sin[#2], 0.0} &, {radii, phis}];
  coordinateAssociation = AssociationThread[vertices -> pts];

  eps = distanceMatrixRadius[pts, kNN, epsilonFactor];
  graph = makeGeometricGraph[pts, eps];

  <|
    "Name" -> name,
    "Graph" -> graph,
    "Points" -> pts,
    "RadiusAssociation" -> radiusAssociation,
    "PhiAssociation" -> phiAssociation,
    "CoordinateAssociation" -> coordinateAssociation,
    "GraphConstruction" -> "matchedFlatStrictRadiusThreshold",
    "GraphK" -> kNN,
    "GraphEpsilon" -> eps,
    "EpsilonFactor" -> epsilonFactor,
    "Sampling" -> "sameRphiAsCurved"
  |>
];

flatSchwarzschildStrict = makeFlatMatchedDataRadiusThreshold[bhSchwarzschild];

strictControlGraphSummary = graphSummaryExtended /@ {
   bhSchwarzschild,
   flatSchwarzschildOriginal,
   flatSchwarzschildStrict
};

Export[
  FileNameJoin[{extraDir, "strict_control_graph_summary_schwarzschild.csv"}],
  strictControlGraphSummary
];

(* Final control-protocol sensitivity table. *)

sensitivityRG = 3;
sensitivityMaxRefsPerBin = All;

ClearAll[runPairSensitivity];

runPairSensitivity[flatData_Association, protocolLabel_String] := Module[
  {bhRows, flatRows, binnedRows, gapRows},
  Print["Running protocol: ", protocolLabel];

  bhRows = diagnosticsForGraphAtRG[
    bhSchwarzschild,
    nRadialBins,
    sensitivityMaxRefsPerBin,
    sensitivityRG
  ];

  flatRows = diagnosticsForGraphAtRG[
    flatData,
    nRadialBins,
    sensitivityMaxRefsPerBin,
    sensitivityRG
  ];

  binnedRows = numericizeRows[aggregateRowsByGraphAndBin[Join[bhRows, flatRows]]];

  Export[
    FileNameJoin[{extraDir, "sensitivity_binned_rows_" <> protocolLabel <> ".csv"}],
    binnedRows
  ];

  gapRows = gapTableFromBinnedRows[
    binnedRows,
    "SchwarzschildFlamm",
    flatData["Name"],
    diagnosticsForRGSweep,
    sensitivityRG
  ];

  Append[#, "Protocol" -> protocolLabel] & /@ gapRows
];

sensitivityOriginalRows =
  runPairSensitivity[flatSchwarzschildOriginal, "originalMatchedFlatKNN"];

sensitivityStrictRows =
  runPairSensitivity[flatSchwarzschildStrict, "strictRadiusMatchedFlat"];

sensitivityControlTable = Join[sensitivityOriginalRows, sensitivityStrictRows];

Export[
  FileNameJoin[{extraDir, "control_protocol_sensitivity_table_schwarzschild.csv"}],
  sensitivityControlTable
];

Export[
  FileNameJoin[{extraDir, "control_protocol_sensitivity_table_schwarzschild.mx"}],
  sensitivityControlTable
];

sensitivityCLogOnly = Select[sensitivityControlTable, #["Diagnostic"] == "CLog" &];

Export[
  FileNameJoin[{extraDir, "control_protocol_sensitivity_CLog_only.csv"}],
  sensitivityCLogOnly
];

Dataset[sensitivityCLogOnly]
