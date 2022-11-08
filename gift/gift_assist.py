from cao_pylib import sbox_model

cnf="(2'+14)(4'+12)(2'+12)(7'+14)(6'+7)(7'+12)(0'+1)(2'+3)(4'+5)(8+9')(9+10')(10+11')(11+15')(14'+15)(12'+13)(5'+14)(3'+5+12)(0'+14)(1+8'+13)(1'+11)(1+13'+15)(3+12+14')(0+3+5+7+12')(0+2+3'+4+7+12'+14')(11'+12+15)(1'+12+13'+14')(13+14+15')(1'+3'+14+15')(1+12+14+15')(0'+13)"
cnf_inv="(6'+8)(4'+10)(5'+8)(7'+10)(1'+8)(1'+10)(0'+1)(2'+3)(4'+5)(6'+7)(12+13')(13+14')(14+15')(9'+15)(8'+9)(10'+11)(3'+15)(3+9+11')(3+11+12')(1+2+5+6+8')(1+2+4+7+10')(2'+10)(9+10+15')(3'+8+9'+10')(2'+9)(8+11+15')(3+8+10+12')(3'+8'+10+11')"

cnf_linear="(0'+1)(2'+3)(4'+5)(6'+7)(7'+10)(8+9')(9+13')(10'+14)(13+15')(11'+15)(5'+10)(10'+12)(10'+11)(1'+12)(11+12'+14)(3+10+14')(3'+15)(1+8'+11)(3'+10+11'+12')(2'+10)(1+14+15')(0'+10)(0+2+5+7+11'+12'+14')"
cnf_linear_inv="(4'+8)(3'+8)(1'+8)(0'+1)(2'+3)(4'+5)(6'+7)(8'+12)(11+13')(13+14')(14+15')(9'+15)(9+12')(5+8+10')(7'+10+12)(8'+10)(7+11'+12)(7+10+15')(5'+15)(6'+8)(5'+7'+8+9')(1+3+4+6+8')"

sbox_model.gen_ine_from_cnf(cnf_linear,'gift/gift_sbox_ine_linear.txt')
sbox_model.gen_ine_from_cnf(cnf_linear_inv,'gift/gift_sbox_ine_linear_inv.txt')

p=[48, 1, 18, 35, 32, 49, 2, 19, 16, 33, 50, 3, 0, 17, 34, 51, 52, 5, 22, 39, 36, 53, 6, 23, 20, 37, 54, 7, 4, 21, 38, 55, 56, 9, 26, 43, 40, 57, 10, 27, 24, 41, 58, 11, 8, 25, 42, 59, 60, 13, 30, 47, 44, 61, 14, 31, 28, 45, 62, 15, 12, 29, 46, 63]
p_inv=[12, 1, 6, 11, 28, 17, 22, 27, 44, 33, 38, 43, 60, 49, 54, 59, 8, 13, 2, 7, 24, 29, 18, 23, 40, 45, 34, 39, 56, 61, 50, 55, 4, 9, 14, 3, 20, 25, 30, 19, 36, 41, 46, 35, 52, 57, 62, 51, 0, 5, 10, 15, 16, 21, 26, 31, 32, 37, 42, 47, 48, 53, 58, 63]