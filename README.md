# ZMath
## 벡터와 행렬과 관련된 기본적인 함수를 제공하는 템플릿 라이브러리.
## 2019 - 08 - 18
부스트에도 이미 QVM이라는 라이브러리가 제공되는 것을 알았지만, 템플릿에 대해서 공부할겸 내 코딩 스타일에 더 어울리는식으로 작성했다. 

기본적으로 제공되는 연산들은 common/zcommonfunctions.h 에 작성했으며 그 곳에 있는 함수들은 몇개를 제외하고 대부분의 타입들에 대해 모두 적용될 수 있게끔 했다.

다만 기본적으로 그 헤더에 제공되는 함수들에 대해서 반환타입을 반드시 명시하도록 작성했다.
본래 의도는 타입을 지정하는 방법과 타입을 지정하지 않고 float + int 따위의 연산에서 float로 반환타입을 정해주는 방식을 택하고 싶었다.
하지만 템플릿의 타입추론 방식때문에 그게 가능하지가 않았다. 설령 그것을 달성한다고 하더라도 직관성이나 성능을 잃을 것 같았다.
그래서 그 소소한 불편함은 감수하는 식으로 하도록 했다. 
특수화된 클래스들 vec2, vec3 등에서는 반환타입을 명시하지않게 다시 특수화를 해두어 불편함을 조금이나마 덜었다.
  ex) int a = clamp<int>(2,0,1);   vec2 v = clamp<float>(vec2{0,1}, 0, 1); //common 헤더에 있는 함수 사용시
  ex) vec2 v = Math::dot<float>(vec2{1,2}, vec2{3,4}); //common 헤더에 있는 함수 사용시
  ex) vec2 v = vec2::common::dot(vec2{1,2}, vec2{3,4}); //특수화된 함수 사용시
  ex) auto lenSq = v.dot(v); //멤버 함수 사용시
  
이 라이브러리에서 제공하는 기능.

## 1. 편하게 사용할 수 있는 이니셜라이져.
'''  vec2 v = {1,2};
  mat2 m = {v, v};
  m = {1, 2, 3, 4};
'''
## 2. 연산자 오버로딩
  vec2 v = v1 + v2;
  mat2 m = {v, v};
  v *= m;
  bool b = v == m[0];

## 3. 상수 표현식의 사용 가능.
  constexpr vec2 unit_x = vec2{1, 0};
  constexpr vec_base<float, 10> = vec_base<float, 10>{1,2,3,4,5,6,7,8,9};

## 4. 잘못된 표현식을 사용하는 경우.
동적 어설션을 정적 어설션으로 바꿀 수 있게 매크로로 처리하였는데 정적 어설션이 가끔 되려 어디서 틀렸는지 찾기가 어려움을 느꼈다.
그래서, 에러를 찾아내기 최대한 편하게 동적 어설션을 하도록 해두었다.
  vec2 v;
  mat3 m;
  v *= m; //기본 값은 동적 어설션.
  
## 5. 나눗셈에서..
나눗셈은 모두 div<type>(param0, param1); 함수를 사용했다.
division by zero에 대한 핸들링을 이곳에서 추가 할 수 있도록 하는 의도다.
추가는 하지 않았지만 여지는 남겨두었다.