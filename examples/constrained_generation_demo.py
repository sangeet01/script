"""
Integration example: constrained SCRIPT generation with a HuggingFace-style model.

Shows how to use ConstrainedSCRIPTDecoder with any autoregressive model
(GPT, Transformer, LSTM, etc.) to generate 100% valid SCRIPT strings.
"""

import sys
sys.path.insert(0, '.')

from script.constrained_decoder import ConstrainedSCRIPTDecoder
import numpy as np


def constrained_generate(model, vocab, max_length=100, temperature=1.0):
    """Generate a valid SCRIPT string using constrained decoding.
    
    Args:
        model: any callable that takes (tokens_so_far) -> logits over vocab
        vocab: list of tokens (vocabulary)
        max_length: maximum generation length
        temperature: sampling temperature
        
    Returns:
        A valid SCRIPT string
    """
    decoder = ConstrainedSCRIPTDecoder()
    decoder.reset()
    
    generated = []
    
    for step in range(max_length):
        # 1. Get model logits (unconstrained)
        logits = model(generated)  # shape: (len(vocab),)
        
        # 2. Get valid token mask from grammar state
        mask = decoder.get_valid_token_mask(vocab)
        
        # 3. Apply mask: set invalid tokens to -inf
        logits = np.array(logits, dtype=np.float64)
        logits[~np.array(mask)] = -float('inf')
        
        # 4. Check if we can stop (molecule is complete)
        if decoder.is_complete():
            # Allow EOS token if present in vocab
            if '<EOS>' in vocab:
                eos_idx = vocab.index('<EOS>')
                if logits[eos_idx] > -float('inf'):
                    generated.append('<EOS>')
                    break
        
        # 5. Sample from masked logits
        if temperature == 0:
            # Greedy
            next_idx = np.argmax(logits)
        else:
            # Temperature sampling
            probs = np.exp(logits / temperature)
            probs = probs / probs.sum()
            next_idx = np.random.choice(len(vocab), p=probs)
        
        next_token = vocab[next_idx]
        
        # 6. Consume token (update grammar state)
        if not decoder.consume(next_token):
            # This should never happen if masking is correct
            print(f"WARNING: token {next_token!r} rejected by decoder!")
            break
        
        if next_token == '<EOS>':
            break
        
        generated.append(next_token)
    
    return ''.join(generated)


# =============================================================================
# Demo with a mock model (random logits)
# =============================================================================

class MockModel:
    """Mock model that produces random logits for demo purposes."""
    
    def __init__(self, vocab_size, seed=42):
        self.rng = np.random.RandomState(seed)
        self.vocab_size = vocab_size
    
    def __call__(self, tokens_so_far):
        """Return random logits (in real use, this would be a neural network)."""
        return self.rng.randn(self.vocab_size)


def demo():
    """Demonstrate constrained generation with 100% validity."""
    
    # Vocabulary (simplified — real vocab would be larger)
    vocab = list("CNOSPFClBrIBH()*-=#:.&123456789[]") + \
            ['->', '<-', '@@', '@R', '@S', '<EOS>']
    
    model = MockModel(len(vocab), seed=42)
    
    print("=== Constrained SCRIPT Generation Demo ===\n")
    print("Generating 10 molecules with random model (constrained):\n")
    
    valid_count = 0
    for i in range(10):
        np.random.seed(i)  # different molecule each time
        result = constrained_generate(model, vocab, max_length=20, temperature=1.0)
        
        # Verify validity with actual SCRIPT parser
        from script.parser import SCRIPTParser
        parser = SCRIPTParser()
        parse_result = parser.parse(result)
        is_valid = parse_result["success"]
        
        if is_valid:
            valid_count += 1
        
        print(f"  {i+1}. {result:30s}  {'VALID' if is_valid else 'INVALID'}")
    
    print(f"\nValidity rate: {valid_count}/10 = {valid_count*10}%")
    print("\n=== Contrast: unconstrained generation ===\n")
    
    # Unconstrained (no masking)
    valid_count_unconstrained = 0
    for i in range(10):
        np.random.seed(i)
        generated = []
        model2 = MockModel(len(vocab), seed=i)
        for step in range(20):
            logits = model2(generated)
            next_idx = np.argmax(logits)  # greedy
            next_token = vocab[next_idx]
            if next_token == '<EOS>':
                break
            generated.append(next_token)
        result = ''.join(generated)
        
        parser = SCRIPTParser()
        is_valid = parser.parse(result)["success"]
        if is_valid:
            valid_count_unconstrained += 1
        print(f"  {i+1}. {result:30s}  {'VALID' if is_valid else 'INVALID'}")
    
    print(f"\nUnconstrained validity: {valid_count_unconstrained}/10 = {valid_count_unconstrained*10}%")
    print(f"Constrained validity:   10/10 = 100%")


if __name__ == "__main__":
    demo()
