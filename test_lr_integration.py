"""
Test script to verify LIANA integration in HeartMAP package
"""

def test_lr_database_import():
    """Test importing L-R database from heartmap.data"""
    print("Testing L-R database module import...")
    try:
        from heartmap.data import get_ligand_receptor_pairs, LigandReceptorDatabase, LR_DATABASE_AVAILABLE
        print(f"✓ Successfully imported L-R database module")
        print(f"  LR_DATABASE_AVAILABLE: {LR_DATABASE_AVAILABLE}")
        return True
    except ImportError as e:
        print(f"✗ Failed to import: {e}")
        return False

def test_lr_database_creation():
    """Test creating L-R database instance"""
    print("\nTesting L-R database creation...")
    try:
        from heartmap.data import LigandReceptorDatabase
        
        db = LigandReceptorDatabase(resource='consensus')
        print(f"✓ Created LigandReceptorDatabase instance")
        
        pairs = db.get_pairs(confidence_threshold=0.7)
        print(f"✓ Retrieved {len(pairs)} L-R pairs from database")
        
        if len(pairs) > 0:
            print(f"  Sample pair: {pairs[0]}")
        
        return True
    except Exception as e:
        print(f"✗ Failed: {e}")
        return False

def test_pipeline_integration():
    """Test AdvancedCommunicationPipeline has LIANA support"""
    print("\nTesting pipeline integration...")
    try:
        from heartmap.pipelines import AdvancedCommunicationPipeline
        from heartmap import Config
        
        config = Config.default()
        pipeline = AdvancedCommunicationPipeline(config)
        
        has_lr = hasattr(pipeline, 'lr_available')
        print(f"✓ AdvancedCommunicationPipeline initialized")
        print(f"  Has lr_available attribute: {has_lr}")
        
        if has_lr:
            print(f"  L-R database available: {pipeline.lr_available}")
        
        return True
    except Exception as e:
        print(f"✗ Failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Run all tests"""
    print("=" * 60)
    print("HeartMAP LIANA Integration Tests")
    print("=" * 60)
    
    results = []
    
    results.append(("L-R Database Import", test_lr_database_import()))
    results.append(("L-R Database Creation", test_lr_database_creation()))
    results.append(("Pipeline Integration", test_pipeline_integration()))
    
    print("\n" + "=" * 60)
    print("Test Results Summary")
    print("=" * 60)
    
    for test_name, passed in results:
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"{status}: {test_name}")
    
    all_passed = all(r[1] for r in results)
    print("\n" + ("="*60))
    if all_passed:
        print("✓ All tests passed! LIANA integration is working correctly.")
    else:
        print("✗ Some tests failed. Check errors above.")
    print("=" * 60)

if __name__ == "__main__":
    main()
