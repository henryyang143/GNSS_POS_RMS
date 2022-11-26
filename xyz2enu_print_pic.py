#!/usr/bin/python
# coding=utf-8



import os
import re
import csv
from math import sin,cos,atan,pi,sqrt
import matplotlib.pyplot as plt
import numpy as np

from matplotlib import font_manager


# 找目录下含该关键词的文件名及其路径
def findfile(filepath,key_word):
    filename = []
    for files in os.listdir(filepath):  # 遍历目录文件名
        if (key_word == 'pos' ):
            if re.match(r'[A-Za-z0-9\.]+\.pos',files):
                filename.append(files)  # 文件名及其路径添加到数组
        elif (key_word == 'snx'):
            if re.match(r'igs[0-9]+\.snx',files):
                filename.append(files)  # 文件名及其路径添加到数组
    return filename  # 返回数组

# 经纬度转化为xyz(参考RTKLIB)
def llh2xyz(lat,lon,h):
    RE_WGS84 = 6378137.0
    FE_WGS84 = 1.0/298.257223563
    lat = lat * pi / 180
    lon = lon * pi / 180
    sinp = sin(lat)
    cosp = cos(lat)
    sinl = sin(lon)
    cosl = cos(lon)
    e2 = FE_WGS84*(2.0-FE_WGS84)
    v = RE_WGS84 / sqrt(1.0 - e2 * sinp * sinp)
    x = (v+h)*cosp*cosl
    y = (v+h)*cosp*sinl
    z = (v*(1.0-e2)+h)*sinp
    return [x,y,z]


# 从pos文件中读取坐标结果，并以xyz坐标形式输出
def readpos(posfilepath) :
    f = open(posfilepath,encoding='gb18030', errors='ignore')
    ln = f.readline()
    while ln:
        if '%  GPST' in ln:
            if 'x-ecef(m)' in ln:
                posmode = 'xyz'
            elif 'latitude(deg)' in ln:
                posmode = 'llh'
            break
        ln = f.readline()
    xyz = []
    while ln:
        ln = f.readline()
        if not ln:
            break
        #content = ln.split(' ')
        if posmode == 'xyz':
            x = float(ln[24:38])
            y = float(ln[39:53])
            z = float(ln[54:69])
            xyz.append([x,y,z])
        elif posmode == 'llh':
            lat = float(ln[24:38])
            lon = float(ln[38:53])
            h = float(ln[53:64])
            xyz.append(llh2xyz(lat,lon,h))
    f.close()
    return xyz


# 根据测站id在snx文件中查找测站精确坐标
def getcrd(siteid,snxfilepath,xyz):
    snxcrd=[]
    if snxfilepath == '':
        print('[WARNING] Not find the snxfile, use the average pos as substitute')
    else:
        f = open(snxfilepath,encoding='gb18030', errors='ignore')
        ln = f.readline()
        while ln:
            ln = f.readline()
            if not ln:
                print('[WARNING] Not find the siteid',siteid,', use the average pos as substitute')
                break
            if ln[14:18]==siteid:
                snxcrd.append(float(ln[47:68]))
                ln = f.readline()
                snxcrd.append(float(ln[47:68]))
                ln = f.readline()
                snxcrd.append(float(ln[47:68]))
                break
        f.close()
    if snxcrd == []:
        x = 0
        y = 0
        z = 0
        for i in range(len(xyz)):
            x = x + xyz[i][0]
            y = y + xyz[i][1]
            z = z + xyz[i][2]
        snxcrd = [x/len(xyz),y/len(xyz),z/len(xyz)]
    else:
        print('[INFO] Find the sitecrd in the snx file')
    print('[INFO] The sitecrd is', snxcrd)
    return snxcrd


# xyz转换为llh（经纬度）
def xyz2llh(ecef):
    aell = 6378137 # WGS - 84 not ITRF
    fell = 1 / 298.257223563
    deg = pi / 180
    u = ecef[0]
    v = ecef[1]
    w = ecef[2]
    esq = 2*fell-fell*fell
    if w == 0:
        lat = 0
    else :
        lat0 = atan(w/(1-esq)*sqrt(u*u+v*v))
        j = 0
        delta = 10^6
        limit = 0.000001/3600*deg
        while delta > limit:
            N = aell / sqrt(1 - esq * sin(lat0)*sin(lat0))
            lat = atan((w / sqrt(u*u + v*v)) * (1 + (esq * N * sin(lat0) / w)))
            delta = abs(lat0 - lat)
            lat0 = lat
            j = j + 1
            if j > 10:
                break
    if u == 0 and v == 0:
        long = 0
    elif u > 0 and v >= 0:
        long = atan(v/u)
    elif u < 0:
        long = pi + atan(v/u)
    elif u > 0 and v < 0:
        long = 2*pi
    elif u == 0 and v > 0:
        long = pi/2
    elif u == 0 and v < 0:
        long = 3*pi/2
    h=(sqrt(u*u+v*v)/cos(lat))-N
    llh=[lat,long,h]
    return llh


# xyz转换至以测站精确坐标为基准的enu坐标
def xyz2enu(xyz,snxcrd):
    enu=[]
    llhcrd = xyz2llh(snxcrd)
    phi = llhcrd[0]
    lam = llhcrd[1]
    sinphi = sin(phi)
    cosphi = cos(phi)
    sinlam = sin(lam)
    coslam = cos(lam)
    for i in range(len(xyz)):
        difxyz=[xyz[i][0]-snxcrd[0],xyz[i][1]-snxcrd[1],xyz[i][2]-snxcrd[2]]
        e = -sinlam*difxyz[0]+coslam*difxyz[1]
        n = -sinphi*coslam*difxyz[0]-sinphi*sinlam*difxyz[1]+cosphi*difxyz[2]
        u =  cosphi*coslam*difxyz[0]+cosphi*sinlam*difxyz[1]+sinphi*difxyz[2]
        enu.append([e,n,u])
    return enu

# 保存enu误差结果为CSV文件
def saveenu(enu,posfilepath):
    f = open(posfilepath[:-4] + '.csv', 'a', newline="")
    csv_write = csv.writer(f)
    for i in range(len(enu)):
        csv_write.writerow([enu[i][0],enu[i][1],enu[i][2]])
    f.close()
    return 0


# 绘制enu误差时间序列图
def plotenu(enu,posfilepath):
    e = []
    n = []
    u = []
    for i in range(len(enu)):
        e.append(enu[i][0])
        n.append(enu[i][1])
        u.append(enu[i][2])
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    plt.axis([0, len(enu)+1, -10, 25])#这里可以调节图片的最大最小值
    ax.plot(range(1,len(enu)+1), e)
    ax.plot(range(1,len(enu)+1), n)
    ax.plot(range(1,len(enu)+1), u)
    ax.grid(True, linestyle='--')
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False
    plt.xlabel('以三十秒为一时间间隔的时间戳编号')
    plt.ylabel('ERROR')
    plt.legend(['E','N','U'])
    #plt.show()
    plt.savefig(posfilepath[:-4] + '.png', dpi=1000)
    plt.close()
    return 0

def cal_RMS(enu):
    e=[]
    n=[]
    u=[]

    for i in range(len(enu)):
        e.append(float(enu[i][0]))
        n.append(float(enu[i][1]))
        u.append(float(enu[i][2]))
    #导出e,n,u的误差
    rms_e=(np.sqrt(np.mean(np.square(e))))
    rms_n=(np.sqrt(np.mean(np.square(n))))
    rms_u=(np.sqrt(np.mean(np.square(u))))
    return rms_e,rms_n,rms_u




if __name__ == '__main__':
    site_id=[]
    rms_e=[]
    rms_n=[]
    rms_u=[]


    # 此处修改文件夹路径，将要处理的pos文件和下载的snx坐标文件统一放于此文件夹下
    filepath = 'C:\\data\\all'
    posfilelist = findfile(filepath,'pos')
    snxfilelist = findfile(filepath,'snx')
    if snxfilelist == []:
        snxfilepath = ''
    else:
        snxfilepath = filepath + '\\' + snxfilelist[0]
        print('[INFO] Find the snxfile :',snxfilelist[0])

    for posfile in posfilelist:

        fileindex = posfilelist.index(posfile) + 1
        posfilepath = filepath + '\\' + posfile
        print('[INFO] Process the '+str(fileindex)+'th posfile :',posfile)
        siteid = posfile[0:4].upper()
        site_id.append(siteid)#获取测站的名称
        print('[INFO] The siteid is :',siteid)
        xyz = readpos(posfilepath)
        snxcrd = getcrd(siteid,snxfilepath,xyz)
        enu = xyz2enu(xyz,snxcrd)
        saveenu(enu,posfilepath)
        plotenu(enu,posfilepath)
        a,b,c=cal_RMS(enu)
        rms_e.append(a)
        rms_n.append(b)
        rms_u.append(c)
        print('[INFO] Finish the '+str(fileindex)+'th posfile processing successfully')
    #绘制rms图并保存

    bar_width = 0.2

    bar_1 = list(range(len(site_id)))
    bar_2 = [i + bar_width for i in bar_1]
    bar_3 = [i + bar_width for i in bar_2]

    # 设置图片尺寸与清晰度
    plt.figure(figsize=(20, 8), dpi=80)

    # 导入数据，绘制条形图
    plt.bar(range(len(site_id)), rms_e, width=bar_width, label='RMS_E')
    plt.bar(bar_2, rms_n, width=bar_width, label='RMS_N')
    plt.bar(bar_3, rms_u, width=bar_width, label='RMS_U')

    # 添加标题
    plt.title('各个测站的ENU坐标的RMS值', size=20)
    # 添加xy轴
    plt.xlabel('测站名称')
    plt.ylabel('RMS值')
    # x轴刻度
    plt.xticks(bar_2, site_id, size=15)
    plt.legend()
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False
    # 展示效果图
    plt.savefig("results.png")