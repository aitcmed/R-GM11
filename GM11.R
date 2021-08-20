# 配置中文字体
showtext::showtext_auto(enable = TRUE)

# GM(1,1)主函数
GM11 <- function(x0, t) { 
    # x0为原始序列
    # t为原始数据加预测数据的长度
    #
    # 一次累加生成序列1-AG0序列(累加序列)
    x1 <- cumsum(x0)
    # 初始化数据集b
    # 长度为length(x0)-1的整数部分个numeric类型且值为0
    b <- numeric(length(x0) - 1)
    # n为length(x0)-1长度
    # 因为需要生成MEAN（紧邻均值）生成序列其长度少1
    n <- length(x0) - 1
    # 生成x1的紧邻均值生成序列
    for (i in 1:n) {
        b[i] <- -(x1[i] + x1[i + 1]) / 2
        # 得序列b，即为x1的紧邻均值生成序列
        b
    }
    D <- numeric(length(x0) - 1)
    D[] <- 1
    # 作B矩阵
    B <- cbind(b, D)
    # B矩阵转置
    BT <- t(B) 
    # 求BT*B得逆
    M <- solve(BT %*% B)
    YN <- numeric(length(x0) - 1)
    YN <- x0[2:length(x0)]
    # 模型的最小二乘估计参数列满足alpha尖
    alpha <- M %*% BT %*% YN
    # 将结果变成一列
    alpha2 <- matrix(alpha, ncol = 1)
    # 得到方程的两个系数
    a <- alpha2[1]
    u <- alpha2[2]
    
    # 输出参数估计值及模拟值
    # 
    # 利用最小二乘法求得参数估计值a,u
    cat( "GM(1,1)参数估计值: ", "\n", 
        "发展系数(-a) =", -a, "\n", 
        "灰色作用量(u) =", u, "\n", "\n") 
    
    # t为原始数据加预测数据的长度
    y <- numeric(length(c(1:t)))
    # 第一个数不变
    y[1] <- x1[1]
    # 将a,u的估计值代入时间响应序列函数计算x1拟合序列y
    for (w in 1:(t - 1)) { 
        y[w + 1] <- (x1[1] - u / a) * exp(-a * w) + u / a
    }
    cat("x(1)的模拟值: ", "\n", y, "\n")
    xy <- numeric(length(y))
    xy[1] <- y[1]
    # 运用后减运算还原得模型输入序列x0预测序列
    for (o in 2:t) {
        xy[o] <- y[o] - y[o - 1]
    }
    cat("x(0)的模拟值: ", "\n", xy, "\n", "\n")
    
    # 计算残差e
    e <- numeric(length(x0))
    for (l in 1:length(x0)) {
        # 得残差序列（未取绝对值）
        e[l] <- x0[l] - xy[l]
    }
    cat("绝对残差: ", "\n", e, "\n")
    # 计算相对误差
    e2 <- numeric(length(x0))
    for (s in 1:length(x0)) {
        # 得相对误差
        e2[s] <- (abs(e[s]) / x0[s])
    }
    cat("相对残差: ", "\n", e2, "\n", "\n")
    cat("残差平方和 =", sum(e^2), "\n")
    cat("平均相对误差 =", 
        sum(e2) / (length(e2) - 1) * 100, "%", "\n")
    cat("相对精度 =", 
        (1 - (sum(e2) / (length(e2) - 1))) * 100, "%", "\n", "\n")
    
    # 后验差比值检验
    avge <- mean(abs(e))
    esum <- sum((abs(e) - avge)^2)
    evar <- esum / (length(e) - 1)
    # 计算残差的均方差se
    se <- sqrt(evar) 
    avgx0 <- mean(x0)
    x0sum <- sum((x0 - avgx0)^2)
    x0var <- x0sum / (length(x0))
    # 计算原序列x0的方差sx
    sx <- sqrt(x0var) 
    # 得验差比值（方差比）
    cv <- se / sx
    # 对后验差比值进行检验，与一般标准进行比较判断预测结果好坏。
    cat("后验差比值(C) =", cv, "\n")
    # 计算小残差概率
    P <- sum((abs(e) - avge) < 0.6745 * sx) / length(e)
    cat("小误差概率(P) =", P, "\n", "\n")
    
    # 判断预测精度等级
    cat("预测精度等级: ", "\n", "")
    if (cv <= 0.35 && P >= 0.95) {
        cat("C <= 0.35, P >= 0.95", "\n", 
            "一级（好）", "\n", 
            "可以预测分析", "\n", "\n")
    } else if ((cv > 0.35 && cv <= 0.5) && (P >= 0.8 && p < 0.95)) {
        cat("0.35 < C <= 0.5, 0.8 <= P < 0.95", "\n", 
            "二级（合格）", "\n", 
            "可以预测分析", "\n", "\n")
    } else if ((cv > 0.5 && cv <= 0.65) && (P >= 0.7 && p < 0.8)) {
        cat("0.5 < C <= 0.65, 0.7 <= P < 0.8", "\n", 
            "三级（勉强合格）", "\n",
            "需要修正精度", "\n", "\n")
    } else if (C > 0.65 && P < 0.7) {
        cat("C > 0.65, P < 0.7", "\n", 
            "四级（不合格）", "\n",
            "需要更换模型", "\n", "\n")
    } else {
        cat("请检查数据和代码")
    }
    
    # 预测模型适用范围
    cat("预测模型适用范围: ", "\n", "")
    if (-a <= 0.3) {
        cat("-a <= 0.3", "\n",
            "适用于短中长期预测", "\n", "\n")
    } else if (-a > 0.3 && -a <= 0.5) {
        cat("0.3 < -a <= 0.5", "\n",
            "可用于短期预测", "\n", "\n")
    } else if (-a > 0.5 && -a <= 0.8) {
        cat("0.5 < -a <= 0.8", "\n",
            "不适用短期预测", "\n", "\n")
    } else if (-a > 0.8 && -a <= 1) {
        cat("0.8 < -a <= 0.1", "\n",
            "应采用残差修正模型", "\n", "\n")
    } else if (-a > 1) {
        cat("-a > 1", "\n",
            "不适合建立GM(1,1)模型", "\n", "\n")
    } else {
        cat("请检查数据和代码")
    }
    
    # 原始序列与预测序列数据框
    df <- data.frame(
        # 原始值空缺用NA填充
        initial = append(initial, rep(NA, num)),
        # 预测值
        forecast = xy
    )
    cat("预测值与实际值对比表: ", "\n")
    print(df)
    
    # 画出输入序列x0的预测序列及x0的比较图像
    plot(
        xy,
        col  = "blue",
        type = "b",
        pch  = 16,
        xlab = "时序",
        ylab = "值"
    )
    points(x0, col = "red", type = "b", pch = 4)
    legend(
        "topleft",
        c("预测", "实际"),
        # title = "",
        pch   = c(16, 4),
        lty   = l,
        col   = c("blue", "red")
    )
}

# 原始数据
initial <- c(65117, 73419, 81990, 90183, 102058, 109668, 121816, 135578, 142618, 156828)

# 预测原始数据以外的num个值
num <- 1

# 运行预测模型
GM11(x0 = initial, t = length(initial) + num)