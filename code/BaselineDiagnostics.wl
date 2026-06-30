(* ::Package:: *)

(* ::Package:: *)

(* BaselineDiagnostics.wl *)

ClearAll[
  shortestPathCountsFromSource,
  shellSizeDiagnostics,
  meanLogMultiplicityDiagnostic,
  localGraphDiagnostics,
  baselineDiagnosticsForVertex,
  baselineDiagnosticsForGraph
];

(* 
  Count the number of shortest paths from one source vertex to all reachable vertices.

  For an unweighted graph, this uses a BFS-style dynamic programming rule:
  sigma[source] = 1
  sigma[v] = sum sigma[u] over predecessors u with dist[u] = dist[v] - 1
*)

shortestPathCountsFromSource[g_Graph, source_] := Module[
  {verts, dist, reachable, ordered, sigma, preds},

  verts = VertexList[g];

  dist = AssociationThread[
    verts,
    GraphDistance[g, source, #] & /@ verts
  ];

  reachable = Select[verts, NumericQ[dist[#]] && dist[#] < Infinity &];

  ordered = SortBy[reachable, dist[#] &];

  sigma = AssociationThread[verts, ConstantArray[0, Length[verts]]];
  sigma[source] = 1;

  Do[
    If[v =!= source,
      preds = Select[
        AdjacencyList[g, v],
        KeyExistsQ[dist, #] && NumericQ[dist[#]] && dist[#] == dist[v] - 1 &
      ];

      sigma[v] = Total[Lookup[sigma, preds, 0]];
    ],
    {v, ordered}
  ];

  sigma
];


shellSizeDiagnostics[g_Graph, source_] := Module[
  {verts, distances, reachableDistances, shellCounts, logShellCounts},

  verts = DeleteCases[VertexList[g], source];

  distances = GraphDistance[g, source, #] & /@ verts;

  reachableDistances = Select[distances, NumericQ[#] && # < Infinity &];

  shellCounts = Counts[reachableDistances];

  logShellCounts = Log[1 + Values[shellCounts]];

  <|
    "MeanLogShellSize" -> If[Length[logShellCounts] > 0, Mean[logShellCounts], Missing["NoShells"]],
    "StdLogShellSize" -> If[Length[logShellCounts] > 1, StandardDeviation[logShellCounts], Missing["TooFewShells"]],
    "MaxShellSize" -> If[Length[shellCounts] > 0, Max[Values[shellCounts]], Missing["NoShells"]],
    "NumberOfShells" -> Length[shellCounts]
  |>
];


meanLogMultiplicityDiagnostic[g_Graph, source_] := Module[
  {sigma, vals, logVals},

  sigma = shortestPathCountsFromSource[g, source];

  vals = DeleteCases[Values[sigma], 0 | 1];

  logVals = Log[vals];

  <|
    "MeanLogMultiplicity" -> If[Length[logVals] > 0, Mean[logVals], 0],
    "StdLogMultiplicity" -> If[Length[logVals] > 1, StandardDeviation[logVals], 0],
    "MaxMultiplicity" -> If[Length[vals] > 0, Max[vals], 1]
  |>
];


localGraphDiagnostics[g_Graph, source_] := Module[
  {deg, clust},

  deg = VertexDegree[g, source];

  clust = Quiet@Check[
    ClusteringCoefficient[g, source],
    Missing["NotAvailable"]
  ];

  <|
    "Degree" -> deg,
    "ClusteringCoefficient" -> clust
  |>
];


baselineDiagnosticsForVertex[g_Graph, source_, radius_: Missing["NoRadius"]] := Module[
  {shell, mult, local},

  shell = shellSizeDiagnostics[g, source];
  mult = meanLogMultiplicityDiagnostic[g, source];
  local = localGraphDiagnostics[g, source];

  Join[
    <|
      "Vertex" -> source,
      "Radius" -> radius
    |>,
    shell,
    mult,
    local
  ]
];


baselineDiagnosticsForGraph[g_Graph, radiusAssociation_: <||>] := Module[
  {verts},

  verts = VertexList[g];

  Dataset[
    baselineDiagnosticsForVertex[
      g,
      #,
      Lookup[radiusAssociation, #, Missing["NoRadius"]]
    ] & /@ verts
  ]
];
