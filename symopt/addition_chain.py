# -*- coding: utf-8 -*-
"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2020 Caleb Bell <Caleb.Andrew.Bell@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

from __future__ import division


shortest_path_costs = {1:1, 2:2, 3:3, 4:3, 5:4, 6:4, 7:5, 8:4, 9:5, 10:5, 11:6, 12:5, 13:6, 14:6,
                       15:6, 16:5, 17:6, 18:6, 19:7, 20:6, 21:7, 22:7, 23:7, 24:6, 25:7, 26:7, 27:7,
                       28:7, 29:8, 30:7, 31:8, 32:6, 33:7, 34:7, 35:8, 36:7, 37:8, 38:8, 39:8, 40:7, 
                       41:8, 42:8, 43:8, 44:8, 45:8, 46:8, 47:9, 48:7, 49:8, 50:8, 51:8, 52:8, 53:9, 
                       54:8, 55:9, 56:8, 57:9, 58:9, 59:9, 60:8, 61:9, 62:9, 63:9, 64:7, 65:8, 66:8,
                       67:9, 68:8, 69:9, 70:9, 71:10, 72:8, 73:9, 74:9, 75:9, 76:9, 77:9, 78:9, 79:10,
                       80:8, 81:9, 82:9, 83:9, 84:9, 85:9, 86:9, 87:10, 88:9, 89:10, 90:9, 91:10, 92:9,
                       93:10, 94:10, 95:10, 96:8, 97:9, 98:9, 99:9, 100:9, 101:10, 102:9, 103:10, 104:9, 
                       105:10, 106:10, 107:10, 108:9, 109:10, 110:10, 111:10, 112:9, 113:10, 114:10, 115:10,
                       116:10, 117:10, 118:10, 119:10, 120:9, 121:10, 122:10, 123:10, 124:10, 125:10, 126:10,
                       127:11, 128:8, 129:9, 130:9, 131:10, 132:9, 133:10, 134:10, 135:10, 136:9, 137:10, 138:10,
                       139:11, 140:10, 141:11, 142:11, 143:11, 144:9, 145:10, 146:10, 147:10, 148:10, 149:10,
                       150:10, 151:11, 152:10, 153:10, 154:10, 155:11, 156:10, 157:11, 158:11, 159:11, 160:9,
                       161:10, 162:10, 163:10, 164:10, 165:10, 166:10, 167:11, 168:10, 169:11, 170:10, 171:11,
                       172:10, 173:11, 174:11, 175:11, 176:10, 177:11, 178:11, 179:11, 180:10, 181:11, 182:11, 
                       183:11, 184:10, 185:11, 186:11, 187:11, 188:11, 189:11, 190:11, 191:12, 192:9, 193:10, 
                       194:10, 195:10, 196:10, 197:11, 198:10, 199:11, 200:10, 201:11, 202:11, 203:11, 204:10, 
                       205:11, 206:11, 207:11, 208:10, 209:11, 210:11, 211:11, 212:11, 213:11, 214:11, 215:11,
                       216:10, 217:11, 218:11, 219:11, 220:11, 221:11, 222:11, 223:12, 224:10, 225:11, 226:11, 
                       227:11, 228:11, 229:11, 230:11, 231:11, 232:11, 233:11, 234:11, 235:12, 236:11, 237:12,
                       238:11, 239:12, 240:10, 241:11, 242:11, 243:11, 244:11, 245:11, 246:11, 247:12, 248:11, 
                       249:11, 250:11, 251:12, 252:11, 253:12, 254:12, 255:11, 256:9, 257:10, 258:10, 259:11, 
                       260:10, 261:11, 262:11, 263:12, 264:10, 265:11, 266:11, 267:12, 268:11, 269:12, 270:11,
                       271:12, 272:10, 273:11, 274:11, 275:12, 276:11, 277:12, 278:12, 279:12, 280:11, 281:11, 
                       282:12, 283:12, 284:12, 285:12, 286:12, 287:12, 288:10, 289:11, 290:11, 291:11, 292:11, 
                       293:11, 294:11, 295:12, 296:11, 297:11, 298:11, 299:12, 300:11, 301:12, 302:12, 303:12,
                       304:11, 305:12, 306:11, 307:12, 308:11, 309:12, 310:12, 311:12, 312:11, 313:12, 314:12, 
                       315:12, 316:12, 317:12, 318:12, 319:12, 320:10, 321:11, 322:11, 323:11}

def addition_chain_length(power):
    return shortest_path_costs[power]

def addition_chains_naive_length(powers):
    '''U
    '''
    # Find set of intermediaries needed for each power
    tot = 0
    for p in powers:
        tot += shortest_path_costs[p]
    return tot
