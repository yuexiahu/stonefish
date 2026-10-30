# 02: 阻力公式的理论文档与变更日志

**What to build:** 理论页的 Drag 小节把 Stonefish 的逐面阻力模型写成可对照标准阻力方程的公式，并在变更日志中记录本次行为变化。

**Blocked by:** None (can start immediately)

**Status:** ready-for-agent

- [x] 理论页 Drag 小节给出逐面形式阻力公式与全局 `0.5·ρ·Cd` 缩放的关系
- [x] 说明迎流面求和等于投影面积，以及无水流直线运动时退化为标准阻力方程
- [x] 说明 `quadratic_drag` / `viscous_drag` 系数的含义与量纲（无量纲、与速度无关）
- [x] 变更日志 1.6 段落新增一条修正说明，措辞与既有条目风格一致
- [x] 文档改动不引入新的 Sphinx 语法或明显格式偏差
