## 修改记录
2025.03.27 
创建文档

## 描述
巡逻路径规划算法，输入信息

## 参数示例：
```json
{
    "id": 1,
    "method": "notice-event",
    "content": {
        "arguments": {
            "taskArea": [
                [121, 22],
                [121, 22],
                [121, 22],
                [121, 22],
                [121, 22]
            ],
            "forbidArea": [
                [
                    [121, 22],
                    [121, 22],
                    [121, 22],
                    [121, 22],
                    [121, 22]
                ],
                [
                    [121, 22],
                    [121, 22],
                    [121, 22],
                    [121, 22]
                ],
            ],
            "vesInfo": [
                {
                    "tid": "03AF0101",
                    "vesPos": [121, 22],
                    "bowAng": 1.0
                },
                {
                    "tid": "03AF0102",
                    "vesPos": [121, 22],
                    "bowAng": 1.0
                }
            ],
            "mustPoint": [
                {
                    "center": [121, 22],
                    "radius": 10,
                    "value": 5,
                    "lastVisit": 50
                }
            ],
            "speed": 10
        }
    }
}
```

## 参数说明

| 参数名         | 类型   |       说明      |
|:--------------|:------:|----------------|
| taskArea | array | 任务区域坐标集合，第一个坐标和最后一个坐标相同 |
| forbidArea | array | 禁止区域坐标集合 |
| vesInfo | array | 各艇信息 |
| tid | string | 艇ID |
| vesPos | array | 艇当前位置 |
| bowAng | float | 艇当前艏向 |
| mustPoint | array | 必过点信息 |
| center | array | 必过点经纬度 |
| radius | float | 必过点径，暂定 |
| value | int | 必过点权重 |
| lastVisit | float | 距离上次访问的时间，单位s |
| speed | float | 指定航速，单位kn |



## 描述
巡逻路径规划算法，输出信息

## 参数示例：
```json
{
    "id": 1,
    "method": "notice-event",
    "content":
    {
        "arguments":
        {
            "vesInfo": [
                {
                    "tid": "01AF0101",
                    "path": [
                        {
                            "shape": "LineString",
                            "points": [
                                {
                                    "coord": [121, 22],
                                    "spd": 10
                                },
                                {
                                    "coord": [121, 22],
                                    "spd": 10
                                }
                            ]
                        },
                        {
                            "shape": "Circle",
                            "points": [
                                {
                                    "coord": [121, 22],
                                    "spd": 10
                                }
                            ],
                            "center": [121, 22],
                            "radius": 100,
                            "reverse": true,
                            "times": 2,
                            "spd": 10
                        }
                    ]
                },
                {
                    "tid": "01AF0102",
                    "path": [
                        {
                            "shape": "LineString",
                            "points": [
                                {
                                    "coord": [121, 22],
                                    "spd": 10
                                },
                                {
                                    "coord": [121, 22],
                                    "spd": 10
                                }
                            ]
                        },
                        {
                            "shape": "Circle",
                            "points": [
                                {
                                    "coord": [121, 22],
                                    "spd": 10
                                }
                            ],
                            "center": [121, 22],
                            "radius": 100,
                            "reverse": true,
                            "times": 2,
                            "spd": 10
                        }
                    ]
                }
            ]
        }
    }
}
```

## 参数说明

| 参数名      | 类型    | 说明 |
|:------------|:-------:|-------------------|
| vesInfo | object | 各艇航路信息 |
| tid | string | 艇ID |
| path | array | 各艇规划路径信息 |
| coord | array | 艇进入环绕切点坐标 |
| spd | float | 驶向各航路点速度 |
| center | array | 环绕圆心坐标 |
| radius | float | 环绕半径 |
| reverse | bool | 环绕方向，true：顺时针，false：逆时针 |
| times | int | 环绕圈数 |
| spd | float | 环绕速度，单位kn |