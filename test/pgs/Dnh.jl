D2hs = [Symel("E", [1 0 0; 0 1 0; 0 0 1]),
    Symel("sigmah", [1 0 0; 0 1 0; 0 0 -1]),
    Symel("i", [-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0]),
    Symel("C_2^1", [-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 1.0]),
    Symel("C_2'(1)", [1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0]),
    Symel("C_2''(1)", [-1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 -1.0]),
    Symel("sigmav_1", [1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 1.0]),
    Symel("sigmad_1", [-1.0 0.0 -0.0; 0.0 1.0 0.0; -0.0 0.0 1.0])]

D3hs = [Symel("E", [1 0 0; 0 1 0; 0 0 1]),
    Symel("sigmah", [1 0 0; 0 1 0; 0 0 -1]),
    Symel("C_3^1", [-0.4999999999999998 -0.8660254037844387 0.0; 0.8660254037844387 -0.4999999999999998 0.0; 0.0 0.0 1.0]),
    Symel("C_3^2", [-0.5000000000000003 0.8660254037844384 0.0; -0.8660254037844384 -0.5000000000000003 0.0; 0.0 0.0 1.0]),
    Symel("C_2'(1)", [1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0]),
    Symel("C_2'(2)", [-0.5000000000000002 -0.8660254037844387 0.0; -0.8660254037844387 0.5000000000000009 0.0; 0.0 0.0 -1.0]),
    Symel("C_2'(3)", [-0.4999999999999991 0.8660254037844394 0.0; 0.8660254037844394 0.4999999999999998 0.0; 0.0 0.0 -1.0]),
    Symel("S_3^1", [-0.4999999999999998 -0.8660254037844387 -0.0; 0.8660254037844387 -0.4999999999999998 0.0; 0.0 0.0 -1.0]),
    Symel("S_3^5", [-0.5000000000000003 0.8660254037844384 0.0; -0.8660254037844384 -0.5000000000000003 -0.0; 0.0 0.0 -1.0]),
    Symel("sigmav_2", [-0.5000000000000009 -0.8660254037844387 0.0; -0.8660254037844387 0.5000000000000002 0.0; 0.0 0.0 1.0]),
    Symel("sigmav_3", [-0.4999999999999998 0.8660254037844394 0.0; 0.8660254037844394 0.4999999999999991 -0.0; 0.0 -0.0 1.0]),
    Symel("sigmav_1", [1.0 -1.2212453270876724e-15 0.0; -1.2212453270876724e-15 -1.0 0.0; 0.0 0.0 1.0])]

D4hs = [Symel("E", [1 0 0; 0 1 0; 0 0 1]),
    Symel("sigmah", [1 0 0; 0 1 0; 0 0 -1]),
    Symel("i", [-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0]),
    Symel("C_4^1", [0.0 -1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 1.0]),
    Symel("C_2^1", [-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 1.0]),
    Symel("C_4^3", [0.0 1.0 0.0; -1.0 0.0 0.0; 0.0 0.0 1.0]),
    Symel("C_2'(1)", [1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0]),
    Symel("C_2'(2)", [-1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 -1.0]),
    Symel("C_2''(3)", [0.0 1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 -1.0]),
    Symel("C_2''(4)", [0.0 -1.0 0.0; -1.0 0.0 0.0; 0.0 0.0 -1.0]),
    Symel("S_4^1", [0.0 -1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 -1.0]),
    Symel("S_4^3", [0.0 1.0 0.0; -1.0 0.0 -0.0; 0.0 0.0 -1.0]),
    Symel("sigmav_2", [1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 1.0]),
    Symel("sigmav_1", [1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 1.0]),
    Symel("sigmad_2", [0.0 -1.0 0.0; -1.0 0.0 0.0; 0.0 0.0 1.0]),
    Symel("sigmad_1", [0.0 1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 1.0])]

D5hs = [Symel("E", [1 0 0; 0 1 0; 0 0 1]),
    Symel("sigmah", [1 0 0; 0 1 0; 0 0 -1]),
    Symel("C_5^1", [0.30901699437494745 -0.9510565162951535 0.0; 0.9510565162951535 0.30901699437494745 0.0; 0.0 0.0 1.0]),
    Symel("C_5^2", [-0.8090169943749473 -0.5877852522924731 0.0; 0.5877852522924731 -0.8090169943749473 0.0; 0.0 0.0 1.0]),
    Symel("C_5^3", [-0.8090169943749475 0.587785252292473 0.0; -0.587785252292473 -0.8090169943749475 0.0; 0.0 0.0 1.0]),
    Symel("C_5^4", [0.30901699437494734 0.9510565162951535 0.0; -0.9510565162951535 0.30901699437494734 0.0; 0.0 0.0 1.0]),
    Symel("C_2'(1)", [1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0]),
    Symel("C_2'(2)", [-0.8090169943749473 0.5877852522924734 0.0; 0.5877852522924734 0.8090169943749481 0.0; 0.0 0.0 -1.0]),
    Symel("C_2'(3)", [0.3090169943749479 -0.951056516295154 0.0; -0.951056516295154 -0.3090169943749471 0.0; 0.0 0.0 -1.0]),
    Symel("C_2'(4)", [0.30901699437494834 0.9510565162951539 0.0; 0.9510565162951539 -0.30901699437494745 0.0; 0.0 0.0 -1.0]),
    Symel("C_2'(5)", [-0.8090169943749475 -0.5877852522924731 0.0; -0.5877852522924731 0.8090169943749481 0.0; 0.0 0.0 -1.0]),
    Symel("S_5^1", [0.30901699437494745 -0.9510565162951535 0.0; 0.9510565162951535 0.30901699437494745 0.0; 0.0 0.0 -1.0]),
    Symel("S_5^7", [-0.8090169943749473 -0.5877852522924731 -0.0; 0.5877852522924731 -0.8090169943749473 0.0; 0.0 0.0 -1.0]),
    Symel("S_5^3", [-0.8090169943749475 0.587785252292473 0.0; -0.587785252292473 -0.8090169943749475 -0.0; 0.0 0.0 -1.0]),
    Symel("S_5^9", [0.30901699437494734 0.9510565162951535 0.0; -0.9510565162951535 0.30901699437494734 0.0; 0.0 0.0 -1.0]),
    Symel("sigmav_2", [-0.8090169943749481 0.5877852522924734 -0.0; 0.5877852522924734 0.8090169943749473 0.0; -0.0 0.0 1.0]),
    Symel("sigmav_3", [0.3090169943749471 -0.951056516295154 0.0; -0.951056516295154 -0.3090169943749479 0.0; 0.0 0.0 1.0]),
    Symel("sigmav_4", [0.30901699437494745 0.9510565162951539 0.0; 0.9510565162951539 -0.30901699437494834 -0.0; 0.0 -0.0 1.0]),
    Symel("sigmav_5", [-0.8090169943749481 -0.5877852522924731 0.0; -0.5877852522924731 0.8090169943749475 0.0; 0.0 0.0 1.0]),
    Symel("sigmav_1", [1.0 -2.2204460492503136e-16 0.0; -2.2204460492503136e-16 -1.0 0.0; 0.0 0.0 1.0])]

D6hs = [Symel("E", [1 0 0; 0 1 0; 0 0 1]),
    Symel("sigmah", [1 0 0; 0 1 0; 0 0 -1]),
    Symel("i", [-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0]),
    Symel("C_6^1", [0.5000000000000001 -0.8660254037844386 0.0; 0.8660254037844386 0.5000000000000001 0.0; 0.0 0.0 1.0]),
    Symel("C_3^1", [-0.4999999999999998 -0.8660254037844388 0.0; 0.8660254037844388 -0.4999999999999998 0.0; 0.0 0.0 1.0]),
    Symel("C_2^1", [-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 1.0]),
    Symel("C_3^2", [-0.5000000000000006 0.8660254037844385 0.0; -0.8660254037844385 -0.5000000000000006 0.0; 0.0 0.0 1.0]),
    Symel("C_6^5", [0.49999999999999944 0.8660254037844392 0.0; -0.8660254037844392 0.49999999999999944 0.0; 0.0 0.0 1.0]),
    Symel("C_2'(1)", [1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0]),
    Symel("C_2'(2)", [-0.4999999999999998 0.8660254037844388 0.0; 0.8660254037844388 0.4999999999999998 0.0; 0.0 0.0 -1.0]),
    Symel("C_2'(3)", [-0.5000000000000004 -0.8660254037844385 0.0; -0.8660254037844385 0.5000000000000007 0.0; 0.0 0.0 -1.0]),
    Symel("C_2''(1)", [0.5000000000000002 0.8660254037844386 0.0; 0.8660254037844386 -0.5000000000000001 0.0; 0.0 0.0 -1.0]),
    Symel("C_2''(2)", [-1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 -1.0]),
    Symel("C_2''(3)", [0.4999999999999998 -0.8660254037844392 0.0; -0.8660254037844392 -0.49999999999999933 0.0; 0.0 0.0 -1.0]),
    Symel("S_6^1", [0.5000000000000001 -0.8660254037844386 0.0; 0.8660254037844386 0.5000000000000001 0.0; 0.0 0.0 -1.0]),
    Symel("S_3^1", [-0.4999999999999998 -0.8660254037844388 -0.0; 0.8660254037844388 -0.4999999999999998 0.0; 0.0 0.0 -1.0]),
    Symel("S_3^5", [-0.5000000000000006 0.8660254037844385 0.0; -0.8660254037844385 -0.5000000000000006 -0.0; 0.0 0.0 -1.0]),
    Symel("S_6^5", [0.49999999999999944 0.8660254037844392 0.0; -0.8660254037844392 0.49999999999999944 0.0; 0.0 0.0 -1.0]),
    Symel("sigmav_2", [-0.5000000000000009 -0.8660254037844387 0.0; -0.8660254037844387 0.5000000000000002 0.0; 0.0 0.0 1.0]),
    Symel("sigmav_3", [-0.4999999999999998 0.8660254037844394 0.0; 0.8660254037844394 0.4999999999999991 -0.0; 0.0 -0.0 1.0]),
    Symel("sigmav_1", [1.0 -1.2212453270876724e-15 0.0; -1.2212453270876724e-15 -1.0 0.0; 0.0 0.0 1.0]),
    Symel("sigmad_2", [-1.0 4.440892098500626e-16 -0.0; 4.440892098500626e-16 1.0 0.0; -0.0 0.0 1.0]),
    Symel("sigmad_3", [0.49999999999999933 -0.8660254037844392 0.0; -0.8660254037844392 -0.4999999999999998 0.0; 0.0 0.0 1.0]),
    Symel("sigmad_1", [0.5000000000000008 0.8660254037844383 0.0; 0.8660254037844383 -0.5000000000000009 -0.0; 0.0 -0.0 1.0])]

D2hcn = ["Ag", "B1g", "B2g", "B3g", "Au", "B1u", "B2u", "B3u"] 
D2hct = [1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0;
         1.0  1.0 -1.0 -1.0  1.0  1.0 -1.0 -1.0;
         1.0 -1.0  1.0 -1.0  1.0 -1.0  1.0 -1.0;
         1.0 -1.0 -1.0  1.0  1.0 -1.0 -1.0  1.0;
         1.0  1.0  1.0  1.0 -1.0 -1.0 -1.0 -1.0;
         1.0  1.0 -1.0 -1.0 -1.0 -1.0  1.0  1.0;
         1.0 -1.0  1.0 -1.0 -1.0  1.0 -1.0  1.0;
         1.0 -1.0 -1.0  1.0 -1.0  1.0  1.0 -1.0]
D3hcn = ["A1'", "A2'", "E'", "A1''", "A2''", "E''"] 
D3hct = D3dct
D4hcn = ["A1g", "A2g", "B1g", "B2g", "Eg", "A1u", "A2u", "B1u", "B2u", "Eu"]
D4hct = [1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0;
         1.0  1.0  1.0 -1.0 -1.0  1.0  1.0  1.0 -1.0 -1.0;
         1.0 -1.0  1.0  1.0 -1.0  1.0 -1.0  1.0  1.0 -1.0;
         1.0 -1.0  1.0 -1.0  1.0  1.0 -1.0  1.0 -1.0  1.0;
         2.0  0.0 -2.0  0.0  0.0  2.0  0.0 -2.0  0.0  0.0;
         1.0  1.0  1.0  1.0  1.0 -1.0 -1.0 -1.0 -1.0 -1.0;
         1.0  1.0  1.0 -1.0 -1.0 -1.0 -1.0 -1.0  1.0  1.0;
         1.0 -1.0  1.0  1.0 -1.0 -1.0  1.0 -1.0 -1.0  1.0;
         1.0 -1.0  1.0 -1.0  1.0 -1.0  1.0 -1.0  1.0 -1.0;
         2.0  0.0 -2.0  0.0  0.0 -2.0  0.0  2.0  0.0  0.0]
D5hcn = ["A1'", "A2'", "E1'", "E2'", "A1''", "A2''", "E1''", "E2''"]
D5hct = D5dct
D6hcn = ["A1g", "A2g", "B1g", "B2g", "E1g", "E2g", "A1u", "A2u", "B1u", "B2u", "E1u", "E2u"]
D6hct = [1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0;
         1.0  1.0  1.0  1.0 -1.0 -1.0  1.0  1.0  1.0  1.0 -1.0 -1.0;
         1.0 -1.0  1.0 -1.0  1.0 -1.0  1.0 -1.0  1.0 -1.0  1.0 -1.0;
         1.0 -1.0  1.0 -1.0 -1.0  1.0  1.0 -1.0  1.0 -1.0 -1.0  1.0;
         2.0  1.0 -1.0 -2.0  0.0  0.0  2.0  1.0 -1.0 -2.0  0.0  0.0;
         2.0 -1.0 -1.0  2.0  0.0  0.0  2.0 -1.0 -1.0  2.0  0.0  0.0;
         1.0  1.0  1.0  1.0  1.0  1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0;
         1.0  1.0  1.0  1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0  1.0  1.0;
         1.0 -1.0  1.0 -1.0  1.0 -1.0 -1.0  1.0 -1.0  1.0 -1.0  1.0;
         1.0 -1.0  1.0 -1.0 -1.0  1.0 -1.0  1.0 -1.0  1.0  1.0 -1.0;
         2.0  1.0 -1.0 -2.0  0.0  0.0 -2.0 -1.0  1.0  2.0  0.0  0.0;
         2.0 -1.0 -1.0  2.0  0.0  0.0 -2.0  1.0  1.0 -2.0  0.0  0.0]