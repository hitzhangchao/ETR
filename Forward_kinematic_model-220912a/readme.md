# 三管ETR扭转运动学模型
> by Zhang Chao
> 2022/7/7

## 正运动学
### 输入：
sheath级 - theta_11, theta_21, theta_31；(arm i, tube 1)
arm1 - theta_12, L_12; theta_13, L_13;    (arm 1, stage i)
arm2 - theta_22, L_22; theta_23, L_23;
arm3 - theta_32, L_32; theta_33, L_33;

### 输出：
T01 - sheath末端位姿;     
T11 - arm1，stage1末端位姿;
T21 - arm2，stage1末端位姿;
T31 - arm3，stage1末端位姿
T12 - arm1，stage2末端位姿;
T22 - arm2，stage2末端位姿;
T32 - arm3，stage2末端位姿
T13 - arm1，stage3末端位姿;
T23 - arm2，stage3末端位姿;
T33 - arm3，stage3末端位姿；

## ETR的形状。


## 运行说明：
clear;clc;

### 进入当前目录
cd('E:\☆ETR\Code\ETR_torsional_rigid_model\Torsional_rigid_model-220707a\');

### 加载CTR物理参数
run("ETR_params.m");

### ETR多臂扭转刚性模型
run("ETR_FK_Rigid.m");

### sheath级，假设三管无位置偏距的扭转柔性模型
run("ETR_sheath_without_offset.m");

### sheath级，基于Lexuan Wang解法的扭转柔性模型
run("ETR_sheath_compliant.m");

### ETR多臂扭转柔性模型
run('ETR_FK_Compliant.m');
