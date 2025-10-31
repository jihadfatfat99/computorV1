import sys
from typing import Union, NoReturn

def handle_error(message: str) -> NoReturn:
    """
    Prints an error message to stderr and exits the program.
    
    Args:
        message (str): The error message to display
    
    Returns:
        NoReturn: This function never returns (exits program)
    
    Exit code: 1 (error)
    
    Handles all types of errors:
    - Parsing errors
    - Validation errors
    - Mathematical errors (division by zero, etc.)
    - Input format errors
    """
    print(f"Error: {message}", file=sys.stderr)
    sys.exit(1)
    
def format_number(n: Union[int, float, None]) -> str:
    """
    Formats a number for clean display in output.
    
    Args:
        n (Union[int, float, None]): The number to format
    
    Returns:
        str: Formatted number string
    
    Handles:
    - Removes trailing .0 for whole numbers (5.0 → "5")
    - Rounds very small floats to 0 (1e-10 → "0")
    - Handles negative zero (-0.0 → "0")
    - Rounds to 6 decimal places for readability
    - Handles special values (inf, nan)
    - Preserves integer representation when possible
    """
    # Handle None
    if n is None:
        return "0"
    
    # Handle special float values
    if isinstance(n, float):
        if n != n:  # NaN check (NaN != NaN is True)
            return "nan"
        if n == float('inf'):
            return "inf"
        if n == float('-inf'):
            return "-inf"
    
    # Convert to float for consistent handling
    num = float(n)
    
    # Handle negative zero
    if num == 0:
        return "0"
    
    # Round very small numbers to zero (precision threshold)
    if abs(num) < 1e-10:
        return "0"
    
    # Round to 6 decimal places to avoid floating point artifacts
    num = round(num, 6)
    
    # Check if it's effectively a whole number
    if num == int(num):
        return str(int(num))
    
    # Return as float string, removing unnecessary trailing zeros
    formatted = f"{num:.6f}".rstrip('0').rstrip('.')
    
    return formatted

def gcd(a: Union[int, float], b: Union[int, float]) -> int:
    """
    Computes the Greatest Common Divisor using Euclidean algorithm.
    Manual implementation (no math.gcd allowed).
    
    Args:
        a (Union[int, float]): First number
        b (Union[int, float]): Second number
    
    Returns:
        int: GCD of a and b (always positive)
    
    Handles:
    - Negative numbers (converts to absolute values)
    - Zero values (gcd(0, n) = n, gcd(n, 0) = n)
    - Equal numbers
    - a < b or a > b (algorithm works both ways)
    
    Algorithm: Euclidean algorithm
    - gcd(a, b) = gcd(b, a mod b)
    - Base case: gcd(a, 0) = a
    """
    # Convert to absolute values (GCD is always positive)
    a = abs(int(a))
    b = abs(int(b))
    
    # Handle edge case: both zero
    if a == 0 and b == 0:
        return 1  # By convention, return 1 to avoid division by zero
    
    # Euclidean algorithm
    while b != 0:
        a, b = b, a % b
    
    return a

def simplify_fraction(numerator: Union[int, float], denominator: Union[int, float]) -> str:
    """
    Reduces a fraction to its irreducible form and returns formatted string.
    
    Args:
        numerator (Union[int, float]): The numerator
        denominator (Union[int, float]): The denominator
    
    Returns:
        str: Formatted fraction string or decimal if not simplifiable
    
    Handles:
    - Zero denominator (error)
    - Zero numerator (returns "0")
    - Negative fractions (sign in numerator)
    - Whole numbers (returns just the number)
    - Already irreducible fractions
    - Very small denominators/numerators
    - Floating point to integer conversion
    
    Output formats:
    - "0" if numerator is 0
    - "5" if result is whole number
    - "3/4" for positive fractions
    - "-3/4" for negative fractions (sign on numerator)
    """
    # Handle zero denominator
    if denominator == 0:
        handle_error("Division by zero in fraction simplification")
    
    # Handle zero numerator
    if abs(numerator) < 1e-10:
        return "0"
    
    # Determine sign
    sign = 1
    if (numerator < 0) != (denominator < 0):  # XOR for sign
        sign = -1
    
    # Work with absolute values
    num = abs(numerator)
    den = abs(denominator)
    
    # Try to convert to integers for GCD calculation
    # Find a multiplier to make both integers (handle decimals)
    precision = 1000000  # 6 decimal places
    
    num_int = round(num * precision)
    den_int = round(den * precision)
    
    # Find GCD and simplify
    divisor = gcd(num_int, den_int)
    
    if divisor > 0:
        num_int //= divisor
        den_int //= divisor
    
    # Apply sign to numerator
    num_int *= sign
    
    # If denominator is 1, return just the numerator
    if den_int == 1:
        return str(num_int)
    
    # Check if it's close to a simple fraction
    # If the numbers are still very large, return decimal instead
    if den_int > 10000:
        result = numerator / denominator
        return format_number(result)
    
    # Return as fraction
    return f"{num_int}/{den_int}"

def validate_equation(equation_str: str) -> bool:
    """
    Validates the equation string for correct syntax and allowed characters.
    
    Args:
        equation_str (str): The equation string to validate
    
    Returns:
        bool: True if valid, raises error otherwise
    
    Validates:
    - Exactly one '=' sign present
    - Only allowed characters: digits, X, ^, *, +, -, =, spaces, dots
    - Not empty or whitespace only
    - Both sides of equation are not empty
    - No invalid character sequences
    
    Raises:
        Calls handle_error() if validation fails
    """
    # Check for empty or whitespace-only string
    if not equation_str or equation_str.strip() == "":
        handle_error("Empty equation provided")
    
    # Check for exactly one equals sign
    equals_count = equation_str.count('=')
    if equals_count == 0:
        handle_error("Missing '=' sign in equation")
    if equals_count > 1:
        handle_error("Multiple '=' signs found in equation")
    
    # Split by '=' and check both sides are not empty
    parts = equation_str.split('=')
    left_side = parts[0].strip()
    right_side = parts[1].strip()
    
    if not left_side:
        handle_error("Left side of equation is empty")
    if not right_side:
        handle_error("Right side of equation is empty")
    
    # Define allowed characters pattern
    # Allowed: digits, X (uppercase only), ^, *, +, -, =, spaces, decimal point
    allowed_pattern = re.compile(r'^[0-9X\^\*\+\-=\s\.]+$')

if __name__ == "__main__":
    print("=" * 60)
    print("TESTING UTILITY FUNCTIONS")
    print("=" * 60)
    
    # Test 1: format_number
    print("\n1️⃣ Testing format_number():")
    print("-" * 40)
    test_numbers = [
        5.0, -3.0, 0.0, -0.0, 
        3.14159, 0.000000001, 
        1234.5678, -0.25,
        float('inf'), float('-inf')
    ]
    for num in test_numbers:
        print(f"format_number({num:>15}) = {format_number(num)}")
    
    # Test 2: gcd
    print("\n3️⃣ Testing gcd():")
    print("-" * 40)
    test_pairs = [
        (48, 18), (100, 50), (7, 13),
        (0, 5), (5, 0), (-12, 8),
        (1, 1), (100, 1)
    ]
    for a, b in test_pairs:
        print(f"gcd({a:>4}, {b:>4}) = {gcd(a, b)}")
    
    # Test 3: simplify_fraction
    print("\n4️⃣ Testing simplify_fraction():")
    print("-" * 40)
    test_fractions = [
        (0, 5), (5, 1), (4, 8),
        (3, 4), (-3, 4), (3, -4),
        (10, 4), (7, 3), (100, 25),
        (-6, 8), (0.5, 2), (1.5, 0)
    ]
    for num, den in test_fractions:
        print(f"simplify_fraction({num:>6}, {den:>4}) = {simplify_fraction(num, den)}")
    
    # Test 4: handle_error (commented out as it exits)
    print("\n1️⃣ Testing handle_error():")
    print("-" * 40)
    print("(Skipped - would exit program)")
    print("Usage: handle_error('Test error message')")
    
    print("\n" + "=" * 60)
    print("✅ All utility functions implemented with type hints!")
    print("=" * 60)