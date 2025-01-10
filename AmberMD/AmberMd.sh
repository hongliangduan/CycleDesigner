#!/bin/bash

# 默认参数
TOPO_FILE=""
COORD_FILE=""
TEMPERATURE=300
SIM_TIME=100

# 显示帮助信息
function show_help {
    echo "使用方法: $0 [选项]"
    echo
    echo "选项:"
    echo "  -p, --topology           指定拓扑文件"
    echo "  -c, --coordinates        指定坐标文件"
    echo "  -temp, --temperature     指定模拟温度 (默认为 300 K)"
    echo "  -time, --simulation-time 指定生产模拟时间 (单位: ns, 默认为 100 ns)"
    echo "  -h, --help               显示此帮助信息"
    echo
    echo "示例:"
    echo "  ./generate_and_run_simulation.sh -p topology.top -c coordinates.crd -temp 300 -time 100"
    exit 0
}

# 解析命令行参数
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -p|--topology) TOPO_FILE="$2"; shift ;;
        -c|--coordinates) COORD_FILE="$2"; shift ;;
        -temp|--temperature) TEMPERATURE="$2"; shift ;;
        -time|--simulation-time) SIM_TIME="$2"; shift ;;
        -h|--help) show_help ;;
        *) echo "未知参数: $1"; exit 1 ;;
    esac
    shift
done

# 检查是否提供了必须的文件参数
if [[ -z "$TOPO_FILE" || -z "$COORD_FILE" ]]; then
    echo "错误: 请提供拓扑文件 (-p) 和坐标文件 (-c)。"
    exit 1
fi

# 使用 cpptraj 从拓扑文件中提取残基数
RESIDUE_COUNT=$(cpptraj -p $TOPO_FILE <<EOF | grep "Total number of residues" | awk '{print $6}'
EOF
)

# 生成 restraintmask 参数
RESTRAINTMASK=":1-${RESIDUE_COUNT}"

# 生成 min.in 文件
cat <<EOL > min.in
min
 &cntrl
  imin=1,maxcyc=1000,ncyc=500,
  cut=8.0,ntb=1,
  ntc=2,ntf=2,
  ntpr=100,
  ntr=1, restraintmask='$RESTRAINTMASK',
  restraint_wt=2.0
 /
EOL

# 生成 heat.in 文件
cat <<EOL > heat.in
heat
 &cntrl
  imin=0,irest=0,ntx=1,
  nstlim=25000,dt=0.002,
  ntc=2,ntf=2,
  cut=8.0, ntb=1,
  ntpr=500, ntwx=500,
  ntt=3, gamma_ln=2.0,
  tempi=0.0, temp0=$TEMPERATURE, ig=-1,
  ntr=1, restraintmask='$RESTRAINTMASK',
  restraint_wt=2.0,
  nmropt=1
 /
 &wt TYPE='TEMP0', istep1=0, istep2=25000,
  value1=0.1, value2=$TEMPERATURE, /
 &wt TYPE='END' /
EOL

# 生成 density.in 文件
cat <<EOL > density.in
density
 &cntrl
  imin=0,irest=1,ntx=5,
  nstlim=25000,dt=0.002,
  ntc=2,ntf=2,
  cut=8.0, ntb=2, ntp=1, taup=1.0,
  ntpr=500, ntwx=500,
  ntt=3, gamma_ln=2.0,
  temp0=$TEMPERATURE, ig=-1,
  ntr=1, restraintmask='$RESTRAINTMASK',
  restraint_wt=2.0,
 /
EOL

# 生成 equil.in 文件
cat <<EOL > equil.in
equil
 &cntrl
  imin=0,irest=1,ntx=5,
  nstlim=250000,dt=0.002,
  ntc=2,ntf=2,
  cut=8.0, ntb=2, ntp=1, taup=2.0,
  ntpr=1000, ntwx=1000,
  ntt=3, gamma_ln=2.0,
  temp0=$TEMPERATURE, ig=-1,
 /
EOL

# 生成 prod.in 文件
NSTLIM=$(($SIM_TIME * 500000))  # 将模拟时间转换为步数，1 ns = 500000 步
cat <<EOL > prod.in
100 ns NPT 生产模拟
 &cntrl
  imin=0,           ! 不进行能量最小化
  irest=1,          ! 从上一次的模拟继续（需要有重启动文件）
  ntx=5,            ! 从重启动文件读取坐标和速度
  nstlim=$NSTLIM,  ! 模拟步数 （模拟时间=(nstlim * dt)/1000=1ns，可以自行调整）
  dt=0.002,         ! 时间步长为2 fs
  ntc=2,            ! 约束含氢键（SHAKE算法）
  ntf=2,            ! 不计算含氢键的力（因为已约束）
  cut=8.0,          ! 非键相互作用截断半径为8.0 Å
  ntb=2,            ! 使用恒压周期性边界条件（NPT）
  ntp=1,            ! 施加各向同性压力耦合
  taup=2.0,         ! 压力弛豫时间为2.0 ps
  ntpr=1000,        ! 每1000步打印一次能量信息
  ntwx=1000,        ! 每1000步输出坐标到轨迹文件
  ntwr=5000,        ! 每5000步写一次重启动文件
  tempi=$TEMPERATURE,      ! 初始温度为$TEMPERATURE K
  temp0=$TEMPERATURE,      ! 目标温度为$TEMPERATURE K
  ntt=3,            ! 使用Langevin动力学进行温度调控
  gamma_ln=1.0,     ! 碰撞频率为1.0 ps^-1
 /
EOL

# 运行能量最小化
pmemd.cuda -O -i min.in -p $TOPO_FILE -c $COORD_FILE -r min.rst -o min.out -ref com_sol.inpcrd -inf min.info

# 运行加热过程
pmemd.cuda -O -i heat.in -p $TOPO_FILE -c min.rst -r heat.rst -o heat.out -ref min.rst -x heat.nc -inf heat.info

# 运行密度平衡
pmemd.cuda -O -i density.in -p $TOPO_FILE -c heat.rst -r density.rst -o density.out -ref heat.rst -x density.nc -inf density.info

# 运行平衡过程
pmemd.cuda -O -i equil.in -p $TOPO_FILE -c density.rst -r equil.rst -o equil.out -x equil.nc -inf equil.info

# 运行生产模拟
pmemd.cuda -O -i prod.in -p $TOPO_FILE -c equil.rst -r prod.rst -o prod.out -x prod.nc -inf prod.info

echo "模拟完成！所有步骤执行完毕。"
