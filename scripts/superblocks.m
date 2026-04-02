superblockRule[Id] = {Subscript[A, 0, 0] -> Subscript[G, 0, 0], 
   Subscript[A, 1, 1] -> 0, 
   Subscript[A, 2, 2] -> 0, 
   Subscript[A, 1, 0] -> 0, 
   Subscript[A, 2, 1] -> 0, 
   Subscript[A, 2, 0] -> 0, 
   Subscript[A, 3, 0] -> 0, Subscript[A, 3, 1] -> 0, 
   Subscript[A, 3, 2] -> 0, Subscript[A, 3, 3] -> 0};

superblockRule[Bp0040] = {Subscript[A, 0, 0] -> (256 Subscript[G, 6, 0])/18375, 
   Subscript[A, 1, 1] -> (256 Subscript[G, 4, 2])/675, 
   Subscript[A, 2, 2] -> Subscript[G, 2, 0], 
   Subscript[A, 1, 0] -> -((128 Subscript[G, 5, 1])/875), 
   Subscript[A, 2, 1] -> -((4 Subscript[G, 3, 1])/3), 
   Subscript[A, 2, 0] -> (16 Subscript[G, 4, 0])/45, 
   Subscript[A, 3, 0] -> 0, Subscript[A, 3, 1] -> 0, 
   Subscript[A, 3, 2] -> 0, Subscript[A, 3, 3] -> 0};

superblockRule[Bp0020] = {Subscript[A, 0, 0] -> Subscript[G, 3, 2]/4, 
   Subscript[A, 1, 1] -> Subscript[G, 1, 0], 
   Subscript[A, 1, 0] -> -Subscript[G, 2, 1], Subscript[A, 2, 0] -> 0,
    Subscript[A, 2, 1] -> 0, Subscript[A, 2, 2] -> 0, 
   Subscript[A, 3, 0] -> 0, Subscript[A, 3, 1] -> 0, 
   Subscript[A, 3, 2] -> 0, Subscript[A, 3, 3] -> 0};

superblockRule[Btwo0200] = {Subscript[A, 0, 
    0] -> (16 Subscript[G, 4, 0])/735 + (512 Subscript[G, 6, 0])/
     56595 + (1024 Subscript[G, 6, 2])/25725 + (
     5120 Subscript[G, 8, 0])/539539, 
   Subscript[A, 1, 
    1] -> (32 Subscript[G, 4, 0])/135 + (512 Subscript[G, 4, 2])/
     945 + (8192 Subscript[G, 6, 2])/25725, 
   Subscript[A, 2, 2] -> (8 Subscript[G, 4, 0])/9, 
   Subscript[A, 1, 
    0] -> -((12 Subscript[G, 3, 1])/35) - (128 Subscript[G, 5, 1])/
     525 - (2304 Subscript[G, 5, 3])/6125 - (1024 Subscript[G, 7, 1])/
     11319, Subscript[A, 2, 
    1] -> -((8 Subscript[G, 3, 1])/3) - (192 Subscript[G, 5, 1])/175, 
   Subscript[A, 2, 0] -> 
    Subscript[G, 2, 0] + (16 Subscript[G, 4, 0])/63 + (
     64 Subscript[G, 4, 2])/45 + (256 Subscript[G, 6, 0])/1225, 
   Subscript[A, 3, 0] -> 0, Subscript[A, 3, 1] -> 0, 
   Subscript[A, 3, 2] -> 0, Subscript[A, 3, 3] -> 0};

superblockRule[Ap0020, j_] := {Subscript[A, 0, 
    0] -> (3 (-1 + j) j Subscript[G, 4 + j, -2 + j])/(
     8 (-1 + 4 j^2)) + (4 j (1 + j) (3 + j) Subscript[G, 4 + j, j])/(
     7 (-1 + 2 j) (3 + 2 j) (7 + 2 j)) + (
     72 (1 + j) (2 + j) (3 + j) (4 + j) Subscript[G, 4 + j, 2 + j])/(
     35 (1 + 2 j) (3 + 2 j) (7 + 2 j) (9 + 2 j)) + (
     64 (2 + j) (3 + j)^3 (4 + j)^2 (5 + j) Subscript[G, 6 + j, 
      2 + j])/(
     7 (3 + 2 j) (5 + 2 j) (7 + 2 j)^3 (9 + 2 j) (11 + 2 j)) + (
     96 (3 + j)^3 (4 + j)^3 (5 + j)^2 (6 + j)^2 Subscript[G, 8 + j, 
      2 + j])/((5 + 2 j) (7 + 2 j)^3 (9 + 2 j)^3 (11 + 2 j)^2 (13 + 
        2 j)), Subscript[A, 1, 1] -> 
    Subscript[G, 2 + j, j] + (
     16 (1 + j) (2 + j) (3 + j) Subscript[G, 4 + j, j])/(
     3 (3 + 2 j)^2 (7 + 2 j)) + (
     64 (6 + 5 j + j^2)^2 Subscript[G, 4 + j, 2 + j])/(
     9 (21 + 20 j + 4 j^2)^2) + (
     48 (1 + j) (2 + j) (3 + j)^2 (4 + j)^2 Subscript[G, 6 + j, 
      j])/((3 + 2 j) (5 + 2 j)^2 (7 + 2 j)^2 (9 + 2 j)) + (
     256 (2 + j) (3 + j)^4 (4 + j)^2 Subscript[G, 6 + j, 2 + j])/(
     3 (3 + 2 j) (5 + 2 j) (7 + 2 j)^4 (9 + 2 j)) + (
     256 (3 + j)^4 (4 + j)^4 Subscript[G, 6 + j, 
      4 + j])/((5 + 2 j)^2 (7 + 2 j)^4 (9 + 2 j)^2), 
   Subscript[A, 2, 2] -> 16/3 Subscript[G, 4 + j, 2 + j], 
   Subscript[A, 1, 
    0] -> -((j Subscript[G, 3 + j, -1 + j])/(1 + 2 j)) - (
     12 (1 + j) (3 + j) Subscript[G, 3 + j, 1 + j])/(
     5 (1 + 2 j) (7 + 2 j)) - (
     6 j (1 + j) (3 + j)^2 Subscript[G, 
      5 + j, -1 + j])/((1 + 2 j) (105 + 142 j + 60 j^2 + 8 j^3)) - (
     48 (2 + j) (3 + j)^2 (137 + 350 j + 270 j^2 + 80 j^3 + 
        8 j^4) Subscript[G, 5 + j, 1 + j])/(
     5 (1 + 2 j) (3 + 2 j) (5 + 2 j)^2 (7 + 2 j)^2 (9 + 2 j)) - (
     192 (2 + j) (3 + j)^4 (4 + j) Subscript[G, 5 + j, 3 + j])/(
     5 (3 + 2 j) (5 + 2 j)^2 (7 + 2 j)^2 (9 + 2 j)) - (
     96 (2 + j) (3 + j)^3 (4 + j)^2 (5 + j)^2 Subscript[G, 7 + j, 
      1 + j])/((5 + 2 j)^2 (7 + 2 j)^3 (9 + 2 j)^2 (11 + 2 j)) - (
     256 (3 + j)^4 (4 + j)^3 (5 + j)^2 Subscript[G, 7 + j, 
      3 + j])/((5 + 2 j)^2 (11 + 2 j) (63 + 32 j + 4 j^2)^3), 
   Subscript[A, 2, 
    1] -> -4 Subscript[G, 3 + j, 1 + j] - (
     32 (2 + j) (3 + j)^2 Subscript[G, 5 + j, 
      1 + j])/((5 + 2 j)^2 (7 + 2 j)) - (
     64 (3 + j)^4 Subscript[G, 5 + j, 3 + j])/(35 + 24 j + 4 j^2)^2, 
   Subscript[A, 2, 
    0] -> (4 (1 + j) Subscript[G, 4 + j, j])/(3 + 2 j) + (
     32 (2 + j) (3 + j) Subscript[G, 4 + j, 2 + j])/(
     3 (3 + 2 j) (7 + 2 j)) + (
     64 (3 + j)^3 (4 + j)^2 Subscript[G, 6 + j, 
      2 + j])/((7 + 2 j)^3 (45 + 28 j + 4 j^2)), 
   Subscript[A, 3, 0] -> 0, Subscript[A, 3, 1] -> 0, 
   Subscript[A, 3, 2] -> 0, Subscript[A, 3, 3] -> 0};

superblockRule[Atwo0100, j_] := {Subscript[A, 0, 
    0] -> -((j Subscript[G, 3 + j, -1 + j])/(4 (1 + 2 j))) - (
     4 (1 + j) (4 + j) Subscript[G, 3 + j, 1 + j])/(
     7 (1 + 2 j) (9 + 2 j)) - (
     5 (-2 + j) (-1 + j) j Subscript[G, 5 + j, -3 + j])/(
     8 (-3 + 2 j) (-1 + 2 j) (1 + 2 j)) - (
     6 (-1 + j) j (1 + j) (2 + j) Subscript[G, 5 + j, -1 + j])/(
     7 (-3 + 2 j) (1 + 2 j) (3 + 2 j) (5 + 2 j)) - (
     144 (2 + j)^2 (3 + j) (-14 + 
        j (5 + j) (5 + 4 j (5 + j))) Subscript[G, 5 + j, 1 + j])/(
     35 (-1 + 2 j) (1 + 2 j) (3 + 2 j) (5 + 2 j) (7 + 2 j) (9 + 
        2 j) (11 + 2 j)) - (
     64 (1 + j) (2 + j)^2 (3 + j)^2 (4 + j) Subscript[G, 5 + j, 
      3 + j])/(
     7 (1 + 2 j) (3 + 2 j) (5 + 2 j)^2 (7 + 2 j) (9 + 2 j)) - (
     96 (2 + j) (3 + j)^2 (4 + j)^2 (5 + j)^2 (6 + j) Subscript[G, 
      7 + j, 1 + j])/(
     7 (5 + 2 j)^2 (7 + 2 j)^2 (9 + 2 j)^2 (11 + 2 j) (13 + 2 j)) - (
     64 (2 + j)^2 (3 + j)^2 (4 + j)^3 (5 + j)^2 Subscript[G, 7 + j, 
      3 + j])/((3 + 2 j) (5 + 2 j)^2 (7 + 2 j)^2 (9 + 2 j)^3 (11 + 
        2 j)) - (
     160 (2 + j) (3 + j) (4 + j)^3 (5 + j)^2 (6 + j)^2 (7 + 
        j)^2 Subscript[G, 9 + j, 
      1 + j])/((5 + 2 j) (7 + 2 j)^2 (9 + 2 j)^3 (11 + 2 j)^2 (13 + 
        2 j)^2 (15 + 2 j)), 
   Subscript[A, 1, 
    1] -> -((2 j Subscript[G, 3 + j, -1 + j])/(1 + 2 j)) - (
     16 (1 + j) (2 + j) Subscript[G, 3 + j, 1 + j])/(
     3 (1 + 2 j) (5 + 2 j)) - (
     8 j (1 + j) (3 + j) (4 + j) Subscript[G, 
      5 + j, -1 + j])/((1 + 2 j) (3 + 2 j) (5 + 2 j) (9 + 2 j)) - (
     32 (2 + j)^2 (3 + j) (822 + 
        j (5 + j) (445 + 52 j (5 + j))) Subscript[G, 5 + j, 1 + j])/(
     9 (5 + 2 j)^3 (7 + 2 j) (9 + 2 j) (3 + 4 j (2 + j))) - (
     256 (2 + j)^2 (3 + j)^3 (4 + j) Subscript[G, 5 + j, 3 + j])/(
     3 (3 + 2 j) (5 + 2 j)^3 (7 + 2 j) (9 + 2 j)) - (
     80 j (1 + j) (2 + j) (4 + j)^2 (5 + j)^2 Subscript[G, 
      7 + j, -1 + 
       j])/((1 + 2 j) (3 + 2 j) (5 + 2 j) (7 + 2 j) (9 + 2 j)^2 (11 + 
        2 j)) - (
     128 (1 + j) (2 + j)^2 (3 + j) (4 + j)^2 (5 + j)^2 Subscript[G, 
      7 + j, 
      1 + j])/((1 + 2 j) (5 + 2 j)^2 (7 + 2 j)^2 (9 + 2 j)^2 (11 + 
        2 j)) - (
     512 (2 + j)^2 (3 + j)^2 (4 + j)^3 (5 + j)^2 Subscript[G, 7 + j, 
      3 + j])/((3 + 2 j) (5 + 2 j)^2 (7 + 2 j)^2 (9 + 2 j)^3 (11 + 
        2 j)), Subscript[A, 2, 
    2] -> -((32 (2 + j) Subscript[G, 5 + j, 1 + j])/(3 (5 + 2 j))), 
   Subscript[A, 1, 0] -> 
    Subscript[G, 2 + j, j] + (
     3 (-1 + j) j Subscript[G, 4 + j, -2 + j])/(
     2 (-1 + 2 j) (1 + 2 j)) + (
     4 (1 + j) (2 + j)^2 (-75 + 44 j (4 + j)) Subscript[G, 4 + j, 
      j])/(5 (-1 + 2 j) (3 + 2 j)^2 (5 + 2 j) (9 + 2 j)) + (
     48 (2 + j)^2 (137 + 2 j (5 + j) (35 + 4 j (5 + j))) Subscript[G, 
      4 + j, 2 + j])/(
     5 (1 + 2 j) (3 + 2 j)^2 (5 + 2 j) (7 + 2 j) (9 + 2 j)) + (
     10 (-1 + j) j (1 + j) (4 + j)^2 Subscript[G, 
      6 + j, -2 + 
       j])/((-1 + 2 j) (1 + 2 j) (3 + 2 j) (7 + 2 j) (9 + 2 j)) + (
     72 (1 + j) (2 + j) (4 + j)^2 (-3 + 2 j (5 + j)) Subscript[G, 
      6 + j, j])/(
     5 (-1 + 2 j) (3 + 2 j) (5 + 2 j) (7 + 2 j) (9 + 2 j) (11 + 
        2 j)) + (
     64 (2 + j)^2 (3 + j)^3 (4 + j)^2 (145 + 44 j (6 + j)) Subscript[
      G, 6 + j, 2 + j])/(
     5 (1 + 2 j) (3 + 2 j) (5 + 2 j)^2 (7 + 2 j)^3 (9 + 2 j) (11 + 
        2 j)) + (
     256 (2 + j)^2 (3 + j)^2 (4 + j)^4 Subscript[G, 6 + j, 
      4 + j])/((3 + 2 j) (5 + 2 j)^2 (7 + 2 j)^3 (9 + 2 j)^2) + (
     160 (1 + j) (2 + j) (3 + j) (4 + j)^2 (5 + j)^2 (6 + 
        j)^2 Subscript[G, 8 + j, 
      j])/((5 + 2 j) (7 + 2 j)^2 (9 + 2 j)^2 (11 + 2 j) (429 + 
        430 j + 108 j^2 + 8 j^3)) + (
     384 (2 + j)^2 (3 + j) (4 + j)^3 (5 + j)^2 (6 + j)^2 Subscript[G, 
      8 + j, 2 + 
       j])/((3 + 2 j) (5 + 2 j) (7 + 2 j)^2 (9 + 2 j)^3 (11 + 
        2 j)^2 (13 + 2 j)), 
   Subscript[A, 2, 
    1] -> (8 (1 + j) Subscript[G, 4 + j, j])/(3 + 2 j) + (
     32 (2 + j)^2 Subscript[G, 4 + j, 
      2 + j])/((3 + 2 j) (5 + 2 j)) + (
     48 (1 + j) (2 + j) (4 + j)^2 Subscript[G, 6 + j, 
      j])/((3 + 2 j) (5 + 2 j) (7 + 2 j) (9 + 2 j)) + (
     128 (2 + j)^2 (3 + j) (4 + j)^2 Subscript[G, 6 + j, 
      2 + j])/((3 + 2 j) (5 + 2 j) (7 + 2 j)^2 (9 + 2 j)), 
   Subscript[A, 2, 
    0] -> -4 Subscript[G, 3 + j, 1 + j] - (
     6 j (1 + j) Subscript[G, 
      5 + j, -1 + j])/((1 + 2 j) (3 + 2 j)) - (
     64 (2 + j) (3 + j (5 + j)) Subscript[G, 5 + j, 1 + j])/(
     3 (1 + 2 j) (5 + 2 j) (9 + 2 j)) - (
     64 (2 + j)^2 (3 + j)^2 Subscript[G, 5 + j, 
      3 + j])/((3 + 2 j) (5 + 2 j)^2 (7 + 2 j)) - (
     96 (2 + j) (3 + j) (4 + j)^2 (5 + j)^2 Subscript[G, 7 + j, 
      1 + j])/((5 + 2 j) (7 + 2 j)^2 (9 + 2 j)^2 (11 + 2 j)), 
   Subscript[A, 3, 0] -> 0, Subscript[A, 3, 1] -> 0, 
   Subscript[A, 3, 2] -> 0, Subscript[A, 3, 3] -> 0};

superblockRule[Azero0000, \[CapitalDelta]_, j_] := {Subscript[A, 0, 0] -> 
    Subscript[G, \[CapitalDelta], j] + (
     16 (1 + j - \[CapitalDelta]) \[CapitalDelta] (3 + \
\[CapitalDelta]) (j + \[CapitalDelta]) Subscript[G, 
      2 + \[CapitalDelta], j])/(
     7 (j - \[CapitalDelta]) (1 + j + \[CapitalDelta]) (-1 + 
        2 \[CapitalDelta]) (7 + 2 \[CapitalDelta])) + (
     16 (-3 + j) (-2 + j) (-1 + j) j (-5 + j - \[CapitalDelta]) (-3 + 
        j - \[CapitalDelta]) (-1 + j - \[CapitalDelta]) (1 + 
        j - \[CapitalDelta]) Subscript[G, 
      4 + \[CapitalDelta], -4 + 
       j])/((-5 + 2 j) (-3 + 2 j) (-1 + 2 j) (1 + 2 j) (-6 + 
        j - \[CapitalDelta]) (-4 + j - \[CapitalDelta]) (-2 + 
        j - \[CapitalDelta]) (j - \[CapitalDelta])) + (
     64 (-2 + j) (-1 + j) j (1 + j) (-3 + j - \[CapitalDelta]) (-1 + 
        j - \[CapitalDelta]) (1 + 
        j - \[CapitalDelta]) (j + \[CapitalDelta]) Subscript[G, 
      4 + \[CapitalDelta], -2 + j])/(
     7 (-5 + 2 j) (-1 + 2 j) (1 + 2 j) (3 + 2 j) (-4 + 
        j - \[CapitalDelta]) (-2 + 
        j - \[CapitalDelta]) (j - \[CapitalDelta]) (1 + 
        j + \[CapitalDelta])) + (288 (-1 + j - \[CapitalDelta]) (1 + 
          j - \[CapitalDelta]) (j + \[CapitalDelta]) (2 + 
          j + \[CapitalDelta]) (1080 - 2106 j - 1533 j^2 + 1146 j^3 + 
          573 j^4 + 2250 \[CapitalDelta] - 4800 j \[CapitalDelta] - 
          3400 j^2 \[CapitalDelta] + 2800 j^3 \[CapitalDelta] + 
          1400 j^4 \[CapitalDelta] + 1575 \[CapitalDelta]^2 - 
          3560 j \[CapitalDelta]^2 - 2480 j^2 \[CapitalDelta]^2 + 
          2160 j^3 \[CapitalDelta]^2 + 1080 j^4 \[CapitalDelta]^2 + 
          450 \[CapitalDelta]^3 - 1040 j \[CapitalDelta]^3 - 
          720 j^2 \[CapitalDelta]^3 + 640 j^3 \[CapitalDelta]^3 + 
          320 j^4 \[CapitalDelta]^3 + 45 \[CapitalDelta]^4 - 
          104 j \[CapitalDelta]^4 - 72 j^2 \[CapitalDelta]^4 + 
          64 j^3 \[CapitalDelta]^4 + 
          32 j^4 \[CapitalDelta]^4) Subscript[G, 4 + \[CapitalDelta], 
        j])/(35 (-3 + 2 j) (-1 + 2 j) (3 + 2 j) (5 + 2 j) (-2 + 
          j - \[CapitalDelta]) (j - \[CapitalDelta]) (1 + 
          j + \[CapitalDelta]) (3 + j + \[CapitalDelta]) (1 + 
          2 \[CapitalDelta]) (3 + 2 \[CapitalDelta]) (7 + 
          2 \[CapitalDelta]) (9 + 2 \[CapitalDelta])) + (
     64 j (1 + j) (2 + j) (3 + j) (1 + 
        j - \[CapitalDelta]) (j + \[CapitalDelta]) (2 + 
        j + \[CapitalDelta]) (4 + j + \[CapitalDelta]) Subscript[G, 
      4 + \[CapitalDelta], 2 + j])/(
     7 (-1 + 2 j) (1 + 2 j) (3 + 2 j) (7 + 
        2 j) (j - \[CapitalDelta]) (1 + j + \[CapitalDelta]) (3 + 
        j + \[CapitalDelta]) (5 + j + \[CapitalDelta])) + (
     16 (1 + j) (2 + j) (3 + j) (4 + j) (j + \[CapitalDelta]) (2 + 
        j + \[CapitalDelta]) (4 + j + \[CapitalDelta]) (6 + 
        j + \[CapitalDelta]) Subscript[G, 4 + \[CapitalDelta], 
      4 + j])/((1 + 2 j) (3 + 2 j) (5 + 2 j) (7 + 2 j) (1 + 
        j + \[CapitalDelta]) (3 + j + \[CapitalDelta]) (5 + 
        j + \[CapitalDelta]) (7 + 
        j + \[CapitalDelta])) + (256 (-3 + j - \[CapitalDelta]) (-1 + 
          j - \[CapitalDelta]) (1 + 
          j - \[CapitalDelta]) (2 + \[CapitalDelta]) (3 + \
\[CapitalDelta])^2 (4 + \[CapitalDelta])^2 (5 + \[CapitalDelta]) (j + \
\[CapitalDelta]) (2 + j + \[CapitalDelta]) (4 + 
          j + \[CapitalDelta]) Subscript[G, 6 + \[CapitalDelta], 
        j])/(7 (-4 + j - \[CapitalDelta]) (-2 + 
          j - \[CapitalDelta]) (j - \[CapitalDelta]) (1 + 
          j + \[CapitalDelta]) (3 + j + \[CapitalDelta]) (5 + 
          j + \[CapitalDelta]) (3 + 2 \[CapitalDelta]) (5 + 
          2 \[CapitalDelta]) (7 + 2 \[CapitalDelta])^2 (9 + 
          2 \[CapitalDelta]) (11 + 2 \[CapitalDelta])) + (256 (-5 + 
          j - \[CapitalDelta]) (-3 + j - \[CapitalDelta]) (-1 + 
          j - \[CapitalDelta]) (1 + 
          j - \[CapitalDelta]) (3 + \[CapitalDelta])^2 (4 + \
\[CapitalDelta])^2 (5 + \[CapitalDelta])^2 (6 + \[CapitalDelta])^2 (j \
+ \[CapitalDelta]) (2 + j + \[CapitalDelta]) (4 + 
          j + \[CapitalDelta]) (6 + j + \[CapitalDelta]) Subscript[G, 
        8 + \[CapitalDelta], 
        j])/((-6 + j - \[CapitalDelta]) (-4 + 
          j - \[CapitalDelta]) (-2 + 
          j - \[CapitalDelta]) (j - \[CapitalDelta]) (1 + 
          j + \[CapitalDelta]) (3 + j + \[CapitalDelta]) (5 + 
          j + \[CapitalDelta]) (7 + j + \[CapitalDelta]) (5 + 
          2 \[CapitalDelta]) (7 + 2 \[CapitalDelta])^2 (9 + 
          2 \[CapitalDelta])^2 (11 + 2 \[CapitalDelta])^2 (13 + 
          2 \[CapitalDelta])), 
   Subscript[A, 1, 
    1] -> (32 (-1 + j) j (-1 + j - \[CapitalDelta]) (1 + 
        j - \[CapitalDelta]) Subscript[G, 
      2 + \[CapitalDelta], -2 + 
       j])/((-1 + 2 j) (1 + 2 j) (-2 + 
        j - \[CapitalDelta]) (j - \[CapitalDelta])) + (
     64 j (1 + j) (1 + 
        j - \[CapitalDelta]) (j + \[CapitalDelta]) Subscript[G, 
      2 + \[CapitalDelta], j])/(
     3 (-1 + 2 j) (3 + 2 j) (j - \[CapitalDelta]) (1 + 
        j + \[CapitalDelta])) + (
     32 (1 + j) (2 + j) (j + \[CapitalDelta]) (2 + 
        j + \[CapitalDelta]) Subscript[G, 2 + \[CapitalDelta], 
      2 + j])/((1 + 2 j) (3 + 2 j) (1 + j + \[CapitalDelta]) (3 + 
        j + \[CapitalDelta])) + (
     256 (-1 + j) j (-3 + j - \[CapitalDelta]) (-1 + 
        j - \[CapitalDelta]) (1 + 
        j - \[CapitalDelta]) (2 + \[CapitalDelta]) (3 + \
\[CapitalDelta]) (j + \[CapitalDelta]) Subscript[G, 
      4 + \[CapitalDelta], -2 + j])/(
     3 (-1 + 2 j) (1 + 2 j) (-4 + j - \[CapitalDelta]) (-2 + 
        j - \[CapitalDelta]) (j - \[CapitalDelta]) (1 + 
        j + \[CapitalDelta]) (3 + 2 \[CapitalDelta]) (7 + 
        2 \[CapitalDelta])) + (
     512 j (1 + j) (-1 + j - \[CapitalDelta]) (1 + 
        j - \[CapitalDelta]) (2 + \[CapitalDelta]) (3 + \
\[CapitalDelta]) (j + \[CapitalDelta]) (2 + 
        j + \[CapitalDelta]) Subscript[G, 4 + \[CapitalDelta], j])/(
     9 (-1 + 2 j) (3 + 2 j) (-2 + 
        j - \[CapitalDelta]) (j - \[CapitalDelta]) (1 + 
        j + \[CapitalDelta]) (3 + j + \[CapitalDelta]) (3 + 
        2 \[CapitalDelta]) (7 + 2 \[CapitalDelta])) + (
     256 (1 + j) (2 + j) (1 + 
        j - \[CapitalDelta]) (2 + \[CapitalDelta]) (3 + \
\[CapitalDelta]) (j + \[CapitalDelta]) (2 + j + \[CapitalDelta]) (4 + 
        j + \[CapitalDelta]) Subscript[G, 4 + \[CapitalDelta], 
      2 + j])/(
     3 (1 + 2 j) (3 + 2 j) (j - \[CapitalDelta]) (1 + 
        j + \[CapitalDelta]) (3 + j + \[CapitalDelta]) (5 + 
        j + \[CapitalDelta]) (3 + 2 \[CapitalDelta]) (7 + 
        2 \[CapitalDelta])) + (512 (-1 + j) j (-5 + 
          j - \[CapitalDelta]) (-3 + j - \[CapitalDelta]) (-1 + 
          j - \[CapitalDelta]) (1 + 
          j - \[CapitalDelta]) (3 + \[CapitalDelta])^2 (4 + \
\[CapitalDelta])^2 (j + \[CapitalDelta]) (2 + 
          j + \[CapitalDelta]) Subscript[G, 
        6 + \[CapitalDelta], -2 + 
         j])/((-1 + 2 j) (1 + 2 j) (-6 + j - \[CapitalDelta]) (-4 + 
          j - \[CapitalDelta]) (-2 + 
          j - \[CapitalDelta]) (j - \[CapitalDelta]) (1 + 
          j + \[CapitalDelta]) (3 + j + \[CapitalDelta]) (5 + 
          2 \[CapitalDelta]) (7 + 2 \[CapitalDelta])^2 (9 + 
          2 \[CapitalDelta])) + (1024 j (1 + j) (-3 + 
          j - \[CapitalDelta]) (-1 + j - \[CapitalDelta]) (1 + 
          j - \[CapitalDelta]) (3 + \[CapitalDelta])^2 (4 + \
\[CapitalDelta])^2 (j + \[CapitalDelta]) (2 + 
          j + \[CapitalDelta]) (4 + j + \[CapitalDelta]) Subscript[G, 
        6 + \[CapitalDelta], 
        j])/(3 (-1 + 2 j) (3 + 2 j) (-4 + j - \[CapitalDelta]) (-2 + 
          j - \[CapitalDelta]) (j - \[CapitalDelta]) (1 + 
          j + \[CapitalDelta]) (3 + j + \[CapitalDelta]) (5 + 
          j + \[CapitalDelta]) (5 + 2 \[CapitalDelta]) (7 + 
          2 \[CapitalDelta])^2 (9 + 2 \[CapitalDelta])) + (512 (1 + 
          j) (2 + j) (-1 + j - \[CapitalDelta]) (1 + 
          j - \[CapitalDelta]) (3 + \[CapitalDelta])^2 (4 + \
\[CapitalDelta])^2 (j + \[CapitalDelta]) (2 + 
          j + \[CapitalDelta]) (4 + j + \[CapitalDelta]) (6 + 
          j + \[CapitalDelta]) Subscript[G, 6 + \[CapitalDelta], 
        2 + j])/((1 + 2 j) (3 + 2 j) (-2 + 
          j - \[CapitalDelta]) (j - \[CapitalDelta]) (1 + 
          j + \[CapitalDelta]) (3 + j + \[CapitalDelta]) (5 + 
          j + \[CapitalDelta]) (7 + j + \[CapitalDelta]) (5 + 
          2 \[CapitalDelta]) (7 + 2 \[CapitalDelta])^2 (9 + 
          2 \[CapitalDelta])), 
   Subscript[A, 2, 2] -> (
    128 (-1 + j - \[CapitalDelta]) (1 + 
       j - \[CapitalDelta]) (j + \[CapitalDelta]) (2 + 
       j + \[CapitalDelta]) Subscript[G, 4 + \[CapitalDelta], j])/(
    3 (-2 + j - \[CapitalDelta]) (j - \[CapitalDelta]) (1 + 
       j + \[CapitalDelta]) (3 + j + \[CapitalDelta])), 
   Subscript[A, 1, 
    0] -> -((
      8 j (1 + j - \[CapitalDelta]) Subscript[G, 
       1 + \[CapitalDelta], -1 + 
        j])/((1 + 2 j) (j - \[CapitalDelta]))) - (
     8 (1 + j) (j + \[CapitalDelta]) Subscript[G, 1 + \[CapitalDelta],
       1 + j])/((1 + 2 j) (1 + j + \[CapitalDelta])) - (
     32 (-2 + j) (-1 + j) j (-3 + j - \[CapitalDelta]) (-1 + 
        j - \[CapitalDelta]) (1 + j - \[CapitalDelta]) Subscript[G, 
      3 + \[CapitalDelta], -3 + 
       j])/((-3 + 2 j) (-1 + 2 j) (1 + 2 j) (-4 + 
        j - \[CapitalDelta]) (-2 + 
        j - \[CapitalDelta]) (j - \[CapitalDelta])) - (
     96 j (-1 + j - \[CapitalDelta]) (1 + 
        j - \[CapitalDelta]) (j + \[CapitalDelta]) (-34 + 19 j^2 - 
        52 \[CapitalDelta] + 32 j^2 \[CapitalDelta] - 
        13 \[CapitalDelta]^2 + 8 j^2 \[CapitalDelta]^2) Subscript[G, 
      3 + \[CapitalDelta], -1 + j])/(
     5 (-3 + 2 j) (1 + 2 j) (3 + 2 j) (-2 + 
        j - \[CapitalDelta]) (j - \[CapitalDelta]) (1 + 
        j + \[CapitalDelta]) (1 + 2 \[CapitalDelta]) (7 + 
        2 \[CapitalDelta])) - (96 (1 + j) (1 + 
          j - \[CapitalDelta]) (j + \[CapitalDelta]) (2 + 
          j + \[CapitalDelta]) (-15 + 38 j + 19 j^2 - 
          20 \[CapitalDelta] + 64 j \[CapitalDelta] + 
          32 j^2 \[CapitalDelta] - 5 \[CapitalDelta]^2 + 
          16 j \[CapitalDelta]^2 + 8 j^2 \[CapitalDelta]^2) Subscript[
        G, 3 + \[CapitalDelta], 
        1 + j])/(5 (-1 + 2 j) (1 + 2 j) (5 + 
          2 j) (j - \[CapitalDelta]) (1 + j + \[CapitalDelta]) (3 + 
          j + \[CapitalDelta]) (1 + 2 \[CapitalDelta]) (7 + 
          2 \[CapitalDelta])) - (
     32 (1 + j) (2 + j) (3 + j) (j + \[CapitalDelta]) (2 + 
        j + \[CapitalDelta]) (4 + j + \[CapitalDelta]) Subscript[G, 
      3 + \[CapitalDelta], 
      3 + j])/((1 + 2 j) (3 + 2 j) (5 + 2 j) (1 + 
        j + \[CapitalDelta]) (3 + j + \[CapitalDelta]) (5 + 
        j + \[CapitalDelta])) - (128 (-2 + j) (-1 + j) j (-5 + 
          j - \[CapitalDelta]) (-3 + j - \[CapitalDelta]) (-1 + 
          j - \[CapitalDelta]) (1 + 
          j - \[CapitalDelta]) (3 + \[CapitalDelta])^2 (j + \
\[CapitalDelta]) Subscript[G, 
        5 + \[CapitalDelta], -3 + 
         j])/((-3 + 2 j) (-1 + 2 j) (1 + 2 j) (-6 + 
          j - \[CapitalDelta]) (-4 + j - \[CapitalDelta]) (-2 + 
          j - \[CapitalDelta]) (j - \[CapitalDelta]) (1 + 
          j + \[CapitalDelta]) (5 + 2 \[CapitalDelta]) (7 + 
          2 \[CapitalDelta])) - (384 j (-3 + 
          j - \[CapitalDelta]) (-1 + j - \[CapitalDelta]) (1 + 
          j - \[CapitalDelta]) (3 + \[CapitalDelta])^2 (j + \
\[CapitalDelta]) (2 + j + \[CapitalDelta]) (-99 + 59 j^2 - 
          78 \[CapitalDelta] + 48 j^2 \[CapitalDelta] - 
          13 \[CapitalDelta]^2 + 8 j^2 \[CapitalDelta]^2) Subscript[G,
         5 + \[CapitalDelta], -1 + 
         j])/(5 (-3 + 2 j) (1 + 2 j) (3 + 2 j) (-4 + 
          j - \[CapitalDelta]) (-2 + 
          j - \[CapitalDelta]) (j - \[CapitalDelta]) (1 + 
          j + \[CapitalDelta]) (3 + j + \[CapitalDelta]) (3 + 
          2 \[CapitalDelta]) (5 + 2 \[CapitalDelta]) (7 + 
          2 \[CapitalDelta]) (9 + 2 \[CapitalDelta])) - (384 (1 + 
          j) (-1 + j - \[CapitalDelta]) (1 + 
          j - \[CapitalDelta]) (3 + \[CapitalDelta])^2 (j + \
\[CapitalDelta]) (2 + j + \[CapitalDelta]) (4 + 
          j + \[CapitalDelta]) (-40 + 118 j + 59 j^2 - 
          30 \[CapitalDelta] + 96 j \[CapitalDelta] + 
          48 j^2 \[CapitalDelta] - 5 \[CapitalDelta]^2 + 
          16 j \[CapitalDelta]^2 + 8 j^2 \[CapitalDelta]^2) Subscript[
        G, 5 + \[CapitalDelta], 
        1 + j])/(5 (-1 + 2 j) (1 + 2 j) (5 + 2 j) (-2 + 
          j - \[CapitalDelta]) (j - \[CapitalDelta]) (1 + 
          j + \[CapitalDelta]) (3 + j + \[CapitalDelta]) (5 + 
          j + \[CapitalDelta]) (3 + 2 \[CapitalDelta]) (5 + 
          2 \[CapitalDelta]) (7 + 2 \[CapitalDelta]) (9 + 
          2 \[CapitalDelta])) - (
     128 (1 + j) (2 + j) (3 + j) (1 + 
        j - \[CapitalDelta]) (3 + \[CapitalDelta])^2 (j + \
\[CapitalDelta]) (2 + j + \[CapitalDelta]) (4 + 
        j + \[CapitalDelta]) (6 + j + \[CapitalDelta]) Subscript[G, 
      5 + \[CapitalDelta], 
      3 + j])/((1 + 2 j) (3 + 2 j) (5 + 
        2 j) (j - \[CapitalDelta]) (1 + j + \[CapitalDelta]) (3 + 
        j + \[CapitalDelta]) (5 + j + \[CapitalDelta]) (7 + 
        j + \[CapitalDelta]) (5 + 2 \[CapitalDelta]) (7 + 
        2 \[CapitalDelta])) - (512 j (-5 + j - \[CapitalDelta]) (-3 + 
          j - \[CapitalDelta]) (-1 + j - \[CapitalDelta]) (1 + 
          j - \[CapitalDelta]) (3 + \[CapitalDelta])^2 (4 + \
\[CapitalDelta])^2 (5 + \[CapitalDelta])^2 (j + \[CapitalDelta]) (2 + 
          j + \[CapitalDelta]) (4 + j + \[CapitalDelta]) Subscript[G, 
        7 + \[CapitalDelta], -1 + 
         j])/((1 + 2 j) (-6 + j - \[CapitalDelta]) (-4 + 
          j - \[CapitalDelta]) (-2 + 
          j - \[CapitalDelta]) (j - \[CapitalDelta]) (1 + 
          j + \[CapitalDelta]) (3 + j + \[CapitalDelta]) (5 + 
          j + \[CapitalDelta]) (5 + 2 \[CapitalDelta]) (7 + 
          2 \[CapitalDelta])^2 (9 + 2 \[CapitalDelta])^2 (11 + 
          2 \[CapitalDelta])) - (512 (1 + j) (-3 + 
          j - \[CapitalDelta]) (-1 + j - \[CapitalDelta]) (1 + 
          j - \[CapitalDelta]) (3 + \[CapitalDelta])^2 (4 + \
\[CapitalDelta])^2 (5 + \[CapitalDelta])^2 (j + \[CapitalDelta]) (2 + 
          j + \[CapitalDelta]) (4 + j + \[CapitalDelta]) (6 + 
          j + \[CapitalDelta]) Subscript[G, 7 + \[CapitalDelta], 
        1 + j])/((1 + 2 j) (-4 + j - \[CapitalDelta]) (-2 + 
          j - \[CapitalDelta]) (j - \[CapitalDelta]) (1 + 
          j + \[CapitalDelta]) (3 + j + \[CapitalDelta]) (5 + 
          j + \[CapitalDelta]) (7 + j + \[CapitalDelta]) (5 + 
          2 \[CapitalDelta]) (7 + 2 \[CapitalDelta])^2 (9 + 
          2 \[CapitalDelta])^2 (11 + 2 \[CapitalDelta])), 
   Subscript[A, 2, 
    1] -> -((
      64 j (-1 + j - \[CapitalDelta]) (1 + 
         j - \[CapitalDelta]) (j + \[CapitalDelta]) Subscript[G, 
       3 + \[CapitalDelta], -1 + 
        j])/((1 + 2 j) (-2 + 
         j - \[CapitalDelta]) (j - \[CapitalDelta]) (1 + 
         j + \[CapitalDelta]))) - (
     64 (1 + j) (1 + j - \[CapitalDelta]) (j + \[CapitalDelta]) (2 + 
        j + \[CapitalDelta]) Subscript[G, 3 + \[CapitalDelta], 
      1 + j])/((1 + 2 j) (j - \[CapitalDelta]) (1 + 
        j + \[CapitalDelta]) (3 + j + \[CapitalDelta])) - (
     256 j (-3 + j - \[CapitalDelta]) (-1 + j - \[CapitalDelta]) (1 + 
        j - \[CapitalDelta]) (3 + \[CapitalDelta])^2 (j + \
\[CapitalDelta]) (2 + j + \[CapitalDelta]) Subscript[G, 
      5 + \[CapitalDelta], -1 + 
       j])/((1 + 2 j) (-4 + j - \[CapitalDelta]) (-2 + 
        j - \[CapitalDelta]) (j - \[CapitalDelta]) (1 + 
        j + \[CapitalDelta]) (3 + j + \[CapitalDelta]) (5 + 
        2 \[CapitalDelta]) (7 + 2 \[CapitalDelta])) - (
     256 (1 + j) (-1 + j - \[CapitalDelta]) (1 + 
        j - \[CapitalDelta]) (3 + \[CapitalDelta])^2 (j + \
\[CapitalDelta]) (2 + j + \[CapitalDelta]) (4 + 
        j + \[CapitalDelta]) Subscript[G, 5 + \[CapitalDelta], 
      1 + j])/((1 + 2 j) (-2 + 
        j - \[CapitalDelta]) (j - \[CapitalDelta]) (1 + 
        j + \[CapitalDelta]) (3 + j + \[CapitalDelta]) (5 + 
        j + \[CapitalDelta]) (5 + 2 \[CapitalDelta]) (7 + 
        2 \[CapitalDelta])), 
   Subscript[A, 2, 
    0] -> (16 (1 + 
        j - \[CapitalDelta]) (j + \[CapitalDelta]) Subscript[G, 
      2 + \[CapitalDelta], 
      j])/((j - \[CapitalDelta]) (1 + j + \[CapitalDelta])) + (
     64 (-1 + j) j (-3 + j - \[CapitalDelta]) (-1 + 
        j - \[CapitalDelta]) (1 + 
        j - \[CapitalDelta]) (j + \[CapitalDelta]) Subscript[G, 
      4 + \[CapitalDelta], -2 + 
       j])/((-1 + 2 j) (1 + 2 j) (-4 + j - \[CapitalDelta]) (-2 + 
        j - \[CapitalDelta]) (j - \[CapitalDelta]) (1 + 
        j + \[CapitalDelta])) + (128 (-1 + j - \[CapitalDelta]) (1 + 
          j - \[CapitalDelta]) (j + \[CapitalDelta]) (2 + 
          j + \[CapitalDelta]) (-18 + 45 j + 45 j^2 - 
          15 \[CapitalDelta] + 40 j \[CapitalDelta] + 
          40 j^2 \[CapitalDelta] - 3 \[CapitalDelta]^2 + 
          8 j \[CapitalDelta]^2 + 8 j^2 \[CapitalDelta]^2) Subscript[
        G, 4 + \[CapitalDelta], 
        j])/(3 (-1 + 2 j) (3 + 2 j) (-2 + 
          j - \[CapitalDelta]) (j - \[CapitalDelta]) (1 + 
          j + \[CapitalDelta]) (3 + j + \[CapitalDelta]) (3 + 
          2 \[CapitalDelta]) (7 + 2 \[CapitalDelta])) + (
     64 (1 + j) (2 + j) (1 + 
        j - \[CapitalDelta]) (j + \[CapitalDelta]) (2 + 
        j + \[CapitalDelta]) (4 + j + \[CapitalDelta]) Subscript[G, 
      4 + \[CapitalDelta], 
      2 + j])/((1 + 2 j) (3 + 2 j) (j - \[CapitalDelta]) (1 + 
        j + \[CapitalDelta]) (3 + j + \[CapitalDelta]) (5 + 
        j + \[CapitalDelta])) + (
     256 (-3 + j - \[CapitalDelta]) (-1 + j - \[CapitalDelta]) (1 + 
        j - \[CapitalDelta]) (3 + \[CapitalDelta])^2 (4 + \
\[CapitalDelta])^2 (j + \[CapitalDelta]) (2 + 
        j + \[CapitalDelta]) (4 + j + \[CapitalDelta]) Subscript[G, 
      6 + \[CapitalDelta], 
      j])/((-4 + j - \[CapitalDelta]) (-2 + 
        j - \[CapitalDelta]) (j - \[CapitalDelta]) (1 + 
        j + \[CapitalDelta]) (3 + j + \[CapitalDelta]) (5 + 
        j + \[CapitalDelta]) (5 + 2 \[CapitalDelta]) (7 + 
        2 \[CapitalDelta])^2 (9 + 2 \[CapitalDelta])), 
   Subscript[A, 3, 0] -> 0, Subscript[A, 3, 1] -> 0, 
   Subscript[A, 3, 2] -> 0, Subscript[A, 3, 3] -> 0};